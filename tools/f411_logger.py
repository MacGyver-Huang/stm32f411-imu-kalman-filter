"""
記錄 F411 經 UART 輸出的 CSV
用法:
    python f411_logger.py /dev/ttyUSB0 --seconds 30 --out still.csv

輸出兩個檔案:
    still.csv        原始資料(milli 單位,與 MCU 輸出相同)
    still_roll.csv   換算成 t(s), gyro_rate(deg/s), acc_angle(deg), kf_angle_mcu(deg)
                     可直接給 kalman_sim.py 使用: python kalman_sim.py still_roll.csv
"""
import argparse
import csv
import time
from pathlib import Path

import serial

HEADER_KEY = "t_ms,"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("port", help="例如 /dev/ttyUSB0 或 /dev/ttyACM0")
    p.add_argument("--baud", type=int, default=115200, help="與 CubeMX 的 USART1 設定一致")
    p.add_argument("--seconds", type=float, default=30.0)
    p.add_argument("--out", default="log.csv")
    args = p.parse_args()

    out = Path(args.out)
    roll_out = out.with_name(out.stem + "_roll.csv")

    ser = serial.Serial(args.port, args.baud, timeout=1)
    print(f"等待標頭... (若一直沒反應,按一下板子的 RESET)")

    header = None
    rows = []
    bad = 0
    t_end = None

    while True:
        line = ser.readline().decode(errors="ignore").strip()
        if not line:
            continue
        if line.startswith("#"):
            print(line)                     # 開機訊息、錯誤訊息
            continue
        if header is None:
            if line.startswith(HEADER_KEY):
                header = line.split(",")
                t_end = time.time() + args.seconds
                print(f"開始記錄 {args.seconds:.0f} 秒...")
            continue
        parts = line.split(",")
        if len(parts) != len(header):
            bad += 1
            continue
        try:
            rows.append([int(x) for x in parts])
        except ValueError:
            bad += 1
            continue
        if time.time() >= t_end:
            break

    ser.close()

    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

    idx = {name: i for i, name in enumerate(header)}
    t0 = rows[0][idx["t_ms"]]
    with open(roll_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t", "gyro_rate", "acc_angle", "kf_angle_mcu"])
        for r in rows:
            w.writerow([
                (r[idx["t_ms"]] - t0) / 1000.0,
                r[idx["gx_mdps"]] / 1000.0,
                r[idx["acc_roll_mdeg"]] / 1000.0,
                r[idx["kf_roll_mdeg"]] / 1000.0,
            ])

    kf_ns = [r[idx["kf_ns"]] for r in rows]
    print(f"共 {len(rows)} 筆,格式錯誤 {bad} 行,I2C 錯誤累計 {rows[-1][idx['i2c_err']]}")
    print(f"兩次濾波執行時間: 平均 {sum(kf_ns)/len(kf_ns)/1000:.2f} us, 最大 {max(kf_ns)/1000:.2f} us")
    print(f"已輸出 {out} 與 {roll_out}")


if __name__ == "__main__":
    main()
