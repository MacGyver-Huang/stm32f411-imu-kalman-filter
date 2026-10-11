/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : STM32F411 + MPU6500 姿態估測(二狀態卡爾曼濾波)
  *
  * 修正版重點(全部改在 USER CODE 區塊內,CubeMX 重新產生不會被覆蓋):
  *   [1] 開機檢查 WHO_AM_I,讀錯就停住並快閃 LED
  *   [2] 感測器重置、改用陀螺儀 PLL 時脈、設定量程與數位低通濾波
  *   [3] 校準只用成功讀到的樣本來平均
  *   [4] 讀取失敗時,累積 dt 到下一次成功,並計數錯誤
  *   [5] 用第一筆加速度計角度初始化濾波器
  *   [6] 創新量與角度包裹在 ±180 度
  *   [7] 用 DWT 週期計數器取得微秒級 dt,並量測濾波執行時間
  *   [8] UART 輸出 CSV(整數 milli 單位,不需要浮點 printf)
  *   [9] LED 心跳改用獨立時間戳
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"
#include "dma.h"
#include "i2c.h"
#include "tim.h"
#include "usart.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */
typedef struct {
    float angle;      // 估測角度 (deg)
    float bias;       // 估測的陀螺儀零偏 (deg/s)
    float rate;       // 扣除零偏後的角速度 (deg/s)
    float P[2][2];    // 估測誤差協方差
    float Q_angle;    // 角度過程雜訊
    float Q_bias;     // 零偏過程雜訊
    float R_measure;  // 加速度計角度量測雜訊 (deg^2)
} Kalman_t;

typedef struct {
    int16_t ax, ay, az;
    int16_t gx, gy, gz;
} MPU_Raw_t;
/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */
#define MPU_ADDR               (0x68 << 1)
#define MPU_REG_SMPLRT_DIV     0x19
#define MPU_REG_CONFIG         0x1A
#define MPU_REG_GYRO_CONFIG    0x1B
#define MPU_REG_ACCEL_CONFIG   0x1C
#define MPU_REG_ACCEL_CONFIG2  0x1D
#define MPU_REG_ACCEL_XOUT_H   0x3B
#define MPU_REG_PWR_MGMT_1     0x6B
#define MPU_REG_WHO_AM_I       0x75
#define MPU6500_WHO_AM_I_VAL   0x70   // 依 datasheet,MPU6500 回傳 0x70

#define RAD_TO_DEG             57.2957795f
#define GYRO_LSB_PER_DPS       131.0f  // ±250 dps 量程
#define LOOP_PERIOD_MS         20      // 50 Hz
#define CALIB_SAMPLES          200
/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
static char tx_buf[160];

static Kalman_t kalmanRoll;
static Kalman_t kalmanPitch;

static float gyro_x_offset = 0.0f;   // raw LSB
static float gyro_y_offset = 0.0f;
static float acc_roll_offset = 0.0f; // deg
static float acc_pitch_offset = 0.0f;

static uint32_t i2c_errors = 0;
/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* ---------- 小工具 ---------- */
static float wrap180(float a)                       // [6]
{
    while (a > 180.0f)  a -= 360.0f;
    while (a < -180.0f) a += 360.0f;
    return a;
}

static void UART_Print(const char *s)
{
    HAL_UART_Transmit(&huart1, (uint8_t *)s, (uint16_t)strlen(s), 100);
}

static long to_milli(float x)                       // [8] 浮點轉整數 milli 單位
{
    return (long)lroundf(x * 1000.0f);
}

static void DWT_Init(void)                          // [7]
{
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    DWT->CYCCNT = 0;
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
}

/* ---------- 卡爾曼濾波 ---------- */
void Kalman_Init(Kalman_t *k, float angle0)         // [5]
{
    k->angle = angle0;
    k->bias = 0.0f;
    k->rate = 0.0f;
    k->P[0][0] = 0.0f; k->P[0][1] = 0.0f;
    k->P[1][0] = 0.0f; k->P[1][1] = 0.0f;
    k->Q_angle = 0.001f;
    k->Q_bias = 0.00003f;  // 參數實驗後由 0.003 調降,見 README
    k->R_measure = 0.030f;   // 之後改成靜置量測到的變異數
}

float Kalman_GetAngle(Kalman_t *k, float newAngle, float newRate, float dt)
{
    // 步驟 1:狀態預測
    k->rate = newRate - k->bias;
    k->angle += dt * k->rate;

    // 步驟 2:協方差預測 P = F P F^T + Q
    k->P[0][0] += dt * (dt * k->P[1][1] - k->P[0][1] - k->P[1][0] + k->Q_angle);
    k->P[0][1] -= dt * k->P[1][1];
    k->P[1][0] -= dt * k->P[1][1];
    k->P[1][1] += k->Q_bias * dt;

    // 步驟 3:卡爾曼增益
    float S = k->P[0][0] + k->R_measure;
    float K0 = k->P[0][0] / S;
    float K1 = k->P[1][0] / S;

    // 步驟 4:狀態更新(創新量包裹在 ±180,避免翻越 180 度時跳變) [6]
    float y = wrap180(newAngle - k->angle);
    k->angle = wrap180(k->angle + K0 * y);
    k->bias += K1 * y;

    // 步驟 5:協方差更新 P = (I - K H) P
    float P00 = k->P[0][0];
    float P01 = k->P[0][1];
    k->P[0][0] -= K0 * P00;
    k->P[0][1] -= K0 * P01;
    k->P[1][0] -= K1 * P00;
    k->P[1][1] -= K1 * P01;

    return k->angle;
}

/* ---------- MPU6500 ---------- */
static HAL_StatusTypeDef MPU_Write(uint8_t reg, uint8_t val)
{
    return HAL_I2C_Mem_Write(&hi2c1, MPU_ADDR, reg, 1, &val, 1, 100);
}

// 回傳 true 表示初始化成功;who 回傳讀到的 WHO_AM_I 值   [1][2]
static bool MPU_Init(uint8_t *who)
{
    *who = 0;
    if (HAL_I2C_Mem_Read(&hi2c1, MPU_ADDR, MPU_REG_WHO_AM_I, 1, who, 1, 100) != HAL_OK)
        return false;
    if (*who != MPU6500_WHO_AM_I_VAL)
        return false;

    bool ok = true;
    ok &= MPU_Write(MPU_REG_PWR_MGMT_1, 0x80) == HAL_OK;   // 重置
    HAL_Delay(100);
    ok &= MPU_Write(MPU_REG_PWR_MGMT_1, 0x01) == HAL_OK;   // 喚醒,自動選陀螺儀 PLL 時脈
    HAL_Delay(50);
    ok &= MPU_Write(MPU_REG_SMPLRT_DIV, 4) == HAL_OK;       // 1 kHz / (1+4) = 200 Hz 輸出
    ok &= MPU_Write(MPU_REG_CONFIG, 0x04) == HAL_OK;        // 陀螺儀 DLPF 約 20 Hz(對照 datasheet)
    ok &= MPU_Write(MPU_REG_GYRO_CONFIG, 0x00) == HAL_OK;   // ±250 dps
    ok &= MPU_Write(MPU_REG_ACCEL_CONFIG, 0x00) == HAL_OK;  // ±2 g
    ok &= MPU_Write(MPU_REG_ACCEL_CONFIG2, 0x04) == HAL_OK; // 加速度計 DLPF 約 20 Hz(對照 datasheet)
    HAL_Delay(50);
    return ok;
}

static bool MPU_Read(MPU_Raw_t *r)
{
    uint8_t b[14];
    if (HAL_I2C_Mem_Read(&hi2c1, MPU_ADDR, MPU_REG_ACCEL_XOUT_H, 1, b, 14, 20) != HAL_OK)
        return false;
    r->ax = (int16_t)(b[0] << 8 | b[1]);
    r->ay = (int16_t)(b[2] << 8 | b[3]);
    r->az = (int16_t)(b[4] << 8 | b[5]);
    // b[6], b[7] 是溫度
    r->gx = (int16_t)(b[8] << 8 | b[9]);
    r->gy = (int16_t)(b[10] << 8 | b[11]);
    r->gz = (int16_t)(b[12] << 8 | b[13]);
    return true;
}

static void Acc_Angles(const MPU_Raw_t *r, float *roll, float *pitch)
{
    float ax = r->ax, ay = r->ay, az = r->az;
    *roll  = atan2f(ay, az) * RAD_TO_DEG;
    *pitch = atan2f(-ax, sqrtf(ay * ay + az * az)) * RAD_TO_DEG;
}

// 只用成功讀到的樣本平均   [3]
static bool MPU_Calibrate(void)
{
    float roll_sum = 0, pitch_sum = 0, gx_sum = 0, gy_sum = 0;
    int valid = 0;
    MPU_Raw_t r;

    for (int i = 0; i < CALIB_SAMPLES; i++) {
        if (MPU_Read(&r)) {
            float roll, pitch;
            Acc_Angles(&r, &roll, &pitch);
            roll_sum += roll;
            pitch_sum += pitch;
            gx_sum += r.gx;
            gy_sum += r.gy;
            valid++;
        } else {
            i2c_errors++;
        }
        HAL_Delay(5);
    }
    if (valid < CALIB_SAMPLES / 2)
        return false;

    acc_roll_offset  = roll_sum / valid;
    acc_pitch_offset = pitch_sum / valid;
    gyro_x_offset    = gx_sum / valid;
    gyro_y_offset    = gy_sum / valid;
    return true;
}

// 無法繼續時:快閃 LED,每秒印一次原因
static void Fault_Loop(const char *msg)
{
    uint32_t last = 0;
    while (1) {
        HAL_GPIO_TogglePin(GPIOC, GPIO_PIN_13);
        if (HAL_GetTick() - last >= 1000) {
            UART_Print(msg);
            last = HAL_GetTick();
        }
        HAL_Delay(100);
    }
}
/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{

  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_TIM2_Init();
  MX_USART1_UART_Init();
  MX_I2C1_Init();
  /* USER CODE BEGIN 2 */
  DWT_Init();
  HAL_Delay(50);

  uint8_t who = 0;
  if (!MPU_Init(&who)) {
      snprintf(tx_buf, sizeof(tx_buf), "# ERROR: MPU init failed, WHO_AM_I=0x%02X\r\n", who);
      Fault_Loop(tx_buf);
  }

  UART_Print("# calibrating, keep the board still...\r\n");
  if (!MPU_Calibrate()) {
      Fault_Loop("# ERROR: calibration failed (too many I2C errors)\r\n");
  }

  // 用第一筆加速度計角度初始化濾波器   [5]
  MPU_Raw_t raw;
  float roll0 = 0.0f, pitch0 = 0.0f;
  if (MPU_Read(&raw)) {
      Acc_Angles(&raw, &roll0, &pitch0);
      roll0 -= acc_roll_offset;
      pitch0 -= acc_pitch_offset;
  }
  Kalman_Init(&kalmanRoll, roll0);
  Kalman_Init(&kalmanPitch, pitch0);

  snprintf(tx_buf, sizeof(tx_buf), "# WHO_AM_I=0x%02X clk=%lu\r\n", who, (unsigned long)SystemCoreClock);
  UART_Print(tx_buf);
  UART_Print("t_ms,gx_mdps,gy_mdps,acc_roll_mdeg,acc_pitch_mdeg,kf_roll_mdeg,kf_pitch_mdeg,kf_ns,i2c_err\r\n");

  uint32_t last_tick = HAL_GetTick();
  uint32_t last_led  = last_tick;
  uint32_t last_cyc  = DWT->CYCCNT;
  float pending_dt   = 0.0f;
  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
    uint32_t now = HAL_GetTick();

    if (now - last_tick >= LOOP_PERIOD_MS) {
      last_tick = now;

      // dt 用 DWT 週期計數,解析度遠高於 1 ms   [7]
      uint32_t cyc = DWT->CYCCNT;
      pending_dt += (float)(cyc - last_cyc) / (float)SystemCoreClock;
      last_cyc = cyc;

      if (MPU_Read(&raw)) {
        float gx_rate = (raw.gx - gyro_x_offset) / GYRO_LSB_PER_DPS;
        float gy_rate = (raw.gy - gyro_y_offset) / GYRO_LSB_PER_DPS;
        float acc_roll, acc_pitch;
        Acc_Angles(&raw, &acc_roll, &acc_pitch);
        acc_roll  -= acc_roll_offset;
        acc_pitch -= acc_pitch_offset;

        // 量測兩次濾波共花多少時間   [7]
        uint32_t t0 = DWT->CYCCNT;
        float roll  = Kalman_GetAngle(&kalmanRoll,  acc_roll,  gx_rate, pending_dt);
        float pitch = Kalman_GetAngle(&kalmanPitch, acc_pitch, gy_rate, pending_dt);
        uint32_t kf_cycles = DWT->CYCCNT - t0;
        uint32_t kf_ns = (uint32_t)((uint64_t)kf_cycles * 1000000000ULL / SystemCoreClock);

        pending_dt = 0.0f;   // [4] 成功才歸零;失敗時保留,累積到下一次

        int len = snprintf(tx_buf, sizeof(tx_buf), "%lu,%ld,%ld,%ld,%ld,%ld,%ld,%lu,%lu\r\n",
                           (unsigned long)now,
                           to_milli(gx_rate), to_milli(gy_rate),
                           to_milli(acc_roll), to_milli(acc_pitch),
                           to_milli(roll), to_milli(pitch),
                           (unsigned long)kf_ns, (unsigned long)i2c_errors);
        HAL_UART_Transmit(&huart1, (uint8_t *)tx_buf, (uint16_t)len, 50);
      } else {
        i2c_errors++;   // [4]
      }
    }

    // LED 心跳:獨立時間戳,每 500 ms 翻轉一次   [9]
    if (now - last_led >= 500) {
      HAL_GPIO_TogglePin(GPIOC, GPIO_PIN_13);
      last_led = now;
    }
    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE1);

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSE;
  RCC_OscInitStruct.HSEState = RCC_HSE_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;
  RCC_OscInitStruct.PLL.PLLM = 12;
  RCC_OscInitStruct.PLL.PLLN = 96;
  RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV2;
  RCC_OscInitStruct.PLL.PLLQ = 4;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_3) != HAL_OK)
  {
    Error_Handler();
  }
}

/* USER CODE BEGIN 4 */
void HAL_TIM_PeriodElapsedCallback(TIM_HandleTypeDef *htim) {
  if (htim->Instance == TIM2) {
    // 之後改成計時器中斷驅動的固定週期時使用
  }
}
/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}
#ifdef USE_FULL_ASSERT
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
