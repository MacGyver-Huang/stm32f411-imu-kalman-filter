"""
卡爾曼濾波模擬與參數實驗(與 F411 上的 Kalman_GetAngle 逐行對應)

用法:
    python kalman_sim.py                                  # 模擬資料
    python kalman_sim.py shake_roll.csv                   # 實測資料,預設參數(與韌體相同)
    python kalman_sim.py shake_roll.csv --R 0.3 --Qb 0.0003
    python kalman_sim.py shake_roll.csv --sweep           # 多組參數比較表(不畫圖)
    python kalman_sim.py shake_roll.csv --sweep --R-list 0.002,0.03,0.3 --Qb-list 0.003,0.0003

實測 CSV 欄位: t, gyro_rate, acc_angle[, kf_angle_mcu]   (f411_logger.py 產生的 *_roll.csv)
"""
import argparse
import itertools
from pathlib import Path

import numpy as np

# 與韌體 Kalman_Init 相同的預設值
FW_Q_ANGLE, FW_Q_BIAS, FW_R = 0.001, 0.00003, 0.03


# ---------------------------------------------------------------
# 卡爾曼濾波:刻意寫成和 C 版本一模一樣的純量形式
# ---------------------------------------------------------------
class Kalman:
    def __init__(self, Q_angle, Q_bias, R_measure, angle0=0.0):
        self.angle, self.bias = angle0, 0.0
        self.P = [[0.0, 0.0], [0.0, 0.0]]
        self.Q_angle, self.Q_bias, self.R = Q_angle, Q_bias, R_measure
        self.K = [0.0, 0.0]

    def update(self, new_angle, new_rate, dt):
        P = self.P
        rate = new_rate - self.bias                       # 1. 狀態預測
        self.angle += dt * rate
        P[0][0] += dt * (dt * P[1][1] - P[0][1] - P[1][0] + self.Q_angle)   # 2. 協方差預測
        P[0][1] -= dt * P[1][1]
        P[1][0] -= dt * P[1][1]
        P[1][1] += self.Q_bias * dt
        S = P[0][0] + self.R                              # 3. 卡爾曼增益
        K0, K1 = P[0][0] / S, P[1][0] / S
        y = new_angle - self.angle                        # 4. 狀態更新
        self.angle += K0 * y
        self.bias += K1 * y
        P00, P01 = P[0][0], P[0][1]                       # 5. 協方差更新
        P[0][0] -= K0 * P00
        P[0][1] -= K0 * P01
        P[1][0] -= K1 * P00
        P[1][1] -= K1 * P01
        self.K = [K0, K1]
        return self.angle


def run_kalman(acc, gyro, dt, Qa, Qb, R):
    kf = Kalman(Qa, Qb, R, angle0=acc[0])
    n = acc.size
    est, bias, K0, K1 = (np.zeros(n) for _ in range(4))
    for i in range(n):
        est[i] = kf.update(acc[i], gyro[i], dt[i])
        bias[i] = kf.bias
        K0[i], K1[i] = kf.K
    return est, bias, K0, K1


def complementary(acc, gyro, dt, alpha):
    out = np.zeros_like(acc)
    a = acc[0]
    for i in range(acc.size):
        a = alpha * (a + gyro[i] * dt[i]) + (1 - alpha) * acc[i]
        out[i] = a
    return out


# ---------------------------------------------------------------
# 資料來源
# ---------------------------------------------------------------
def make_sim_data(fs=50.0, duration=40.0, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, duration, 1.0 / fs)
    true = np.zeros_like(t)
    m = (t >= 5) & (t < 15)
    true[m] = 30.0 * np.clip((t[m] - 5) / 2, 0, 1)
    m = (t >= 15) & (t < 30)
    true[m] = 30.0 * np.cos(2 * np.pi * 0.5 * (t[m] - 15))
    m = t >= 30
    true[m] = 30.0 * np.cos(2 * np.pi * 0.5 * 15) * np.exp(-(t[m] - 30))
    true_bias = 2.0 + 0.05 * t
    gyro = np.gradient(true, t) + true_bias + rng.normal(0, 0.3, t.size)
    acc = true + rng.normal(0, 1.5, t.size)
    spikes = rng.random(t.size) < 0.01
    acc[spikes] += rng.normal(0, 15, spikes.sum())
    return dict(t=t, gyro=gyro, acc=acc, true=true, true_bias=true_bias, mcu=None)


def load_csv(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return dict(t=d[:, 0], gyro=d[:, 1], acc=d[:, 2], true=None, true_bias=None,
                mcu=d[:, 3] if d.shape[1] >= 4 else None)


# ---------------------------------------------------------------
# 靜止偵測:靜止時加速度計只量到重力,可以當作真值
# ---------------------------------------------------------------
def detect_still(t, gyro, acc, gyro_thr, acc_std_thr, hold):
    """
    回傳 (still, moving_raw)
    moving_raw[i]: 這一筆看起來在動(角速度大,或加速度計角度在跳)
    still[i]:      往前 hold 秒內都沒有在動,才算真正靜止
    """
    dt_med = float(np.median(np.diff(t)))
    win = max(1, int(round(hold / dt_med)))
    # 加速度計角度的滑動標準差,晃動時會變大
    pad = np.pad(acc, (win - 1, 0), mode="edge")
    roll_std = np.array([pad[i:i + win].std() for i in range(acc.size)])
    moving_raw = (np.abs(gyro) > gyro_thr) | (roll_std > acc_std_thr)
    # 往前 win 筆、往後 win 筆內只要有一筆在動,就不算靜止(離線分析可以看未來)
    m = moving_raw.astype(int)
    k = np.ones(win, dtype=int)
    past = np.convolve(m, k)[:acc.size] > 0
    future = np.convolve(m[::-1], k)[:acc.size][::-1] > 0
    return ~(past | future), moving_raw


def recovery_stats(t, est, acc, still, moving_raw, tol):
    """
    對每一段「動完之後的靜止」:
      恢復時間 = 從最後一次在動,到估測角度和加速度計的差距穩定小於 tol 的時間
      殘差     = 恢復之後,估測角度和加速度計的 RMS 差距
    """
    err = est - acc
    rec_times, resid = [], []
    moved_before = np.cumsum(moving_raw) > 0          # 開頭的校準靜置不算
    i, n = 0, still.size
    while i < n:
        if still[i] and moved_before[i] and (i == 0 or not still[i - 1]):
            s = i
            e = s
            while e < n and still[e]:
                e += 1
            motion_end = np.nonzero(moving_raw[:s])[0][-1]
            seg = np.abs(err[s:e])
            bad = np.nonzero(seg >= tol)[0]
            rec = s if bad.size == 0 else s + bad[-1] + 1
            if rec < e:                                   # 這段內有恢復
                rec_times.append(t[rec] - t[motion_end])
                resid.extend(err[rec:e])
            else:                                         # 整段都沒恢復
                rec_times.append(np.nan)
            i = e
        else:
            i += 1
    return rec_times, resid


def metrics(data, est, bias, args):
    t, acc, gyro, true = data["t"], data["acc"], data["gyro"], data["true"]
    first = (t - t[0]) < args.static_sec
    m = {
        "static_std": float(np.std(est[first])) if first.sum() > 2 else float("nan"),
        "bias_max": float(np.max(np.abs(bias))),
    }
    still, moving_raw = detect_still(t, gyro, acc, args.still_gyro, args.still_acc, args.hold)
    rec, resid = recovery_stats(t, est, acc, still, moving_raw, args.tol)
    m["events"] = len(rec)
    m["recovered"] = int(np.sum(~np.isnan(rec))) if rec else 0
    m["rec_mean"] = float(np.nanmean(rec)) if m["recovered"] else float("nan")
    m["rec_max"] = float(np.nanmax(rec)) if m["recovered"] else float("nan")
    m["resid_rms"] = float(np.sqrt(np.mean(np.square(resid)))) if resid else float("nan")
    m["_still"] = still
    if true is not None:
        m["rms_true"] = float(np.sqrt(np.mean((est - true) ** 2)))
    return m


def fmt(x, unit="", width=7, prec=3):
    return f"{'-':>{width}}" if np.isnan(x) else f"{x:>{width}.{prec}f}{unit}"


def parse_list(s):
    return [float(x) for x in s.split(",") if x.strip()]


def main():
    p = argparse.ArgumentParser(description="卡爾曼濾波模擬與參數實驗")
    p.add_argument("csv", nargs="?", help="*_roll.csv;省略則使用模擬資料")
    p.add_argument("--R", type=float, default=FW_R, help=f"量測雜訊 (預設 {FW_R})")
    p.add_argument("--Qa", type=float, default=FW_Q_ANGLE, help=f"角度過程雜訊 (預設 {FW_Q_ANGLE})")
    p.add_argument("--Qb", type=float, default=FW_Q_BIAS, help=f"零偏過程雜訊 (預設 {FW_Q_BIAS})")
    p.add_argument("--alpha", type=float, default=0.98, help="互補濾波係數 (預設 0.98)")
    p.add_argument("--static-sec", type=float, default=4.0, help="開頭靜置秒數 (預設 4)")
    p.add_argument("--sweep", action="store_true", help="多組參數比較,只印表格")
    p.add_argument("--R-list", default="0.002,0.03,0.3,3", help="sweep 用的 R 清單")
    p.add_argument("--Qb-list", default="0.003,0.0003,0.00003", help="sweep 用的 Qb 清單")
    p.add_argument("--no-show", action="store_true", help="只存圖,不開視窗")
    p.add_argument("--still-gyro", type=float, default=3.0, help="靜止判斷:角速度門檻 deg/s (預設 3)")
    p.add_argument("--still-acc", type=float, default=0.5, help="靜止判斷:加速度計角度滑動標準差門檻 deg (預設 0.5)")
    p.add_argument("--hold", type=float, default=0.5, help="靜止判斷:需持續不動的秒數 (預設 0.5)")
    p.add_argument("--tol", type=float, default=0.5, help="恢復判斷:與加速度計差距小於幾度算恢復 (預設 0.5)")
    args = p.parse_args()

    data = load_csv(args.csv) if args.csv else make_sim_data()
    t, acc, gyro = data["t"], data["acc"], data["gyro"]
    dt = np.diff(t, prepend=t[0] - np.median(np.diff(t)))   # 用實際時間差,和 MCU 一樣

    static = (t - t[0]) < args.static_sec
    if static.sum() > 2:
        print(f"靜置段加速度計角度變異數 (實測 R): {np.var(acc[static]):.4f} deg^2")

    # ---------------- sweep 模式 ----------------
    if args.sweep:
        print(f"\n{'R':>7} {'Qb':>8} | {'靜置雜訊':>8} {'bias偏移':>8} | {'恢復平均':>8} {'恢復最大':>8} {'靜止殘差':>8}")
        print(f"{'':>7} {'':>8} | {'deg':>8} {'deg/s':>8} | {'s':>8} {'s':>8} {'deg':>8}")
        print("-" * 72)
        events_info = None
        for R, Qb in itertools.product(parse_list(args.R_list), parse_list(args.Qb_list)):
            est, bias, _, _ = run_kalman(acc, gyro, dt, args.Qa, Qb, R)
            m = metrics(data, est, bias, args)
            events_info = (m["events"], m["recovered"])
            print(f"{R:>7g} {Qb:>8g} | {fmt(m['static_std'], '', 8, 4)} {fmt(m['bias_max'], '', 8)} |"
                  f" {fmt(m['rec_mean'], '', 8, 2)} {fmt(m['rec_max'], '', 8, 2)} {fmt(m['resid_rms'], '', 8)}")
        print(f"\n偵測到 {events_info[0]} 段「動完後的靜止」(恢復與殘差只在這些時段計算)。")
        print("靜置雜訊越小越平滑;bias 偏移越小越不受干擾;恢復時間越短、殘差越小越好。")
        print("恢復最大顯示 '-' 代表沒有任何一段在靜止期間內恢復,可加大 --tol 或檢查資料。")
        return

    # ---------------- 單組參數 + 畫圖 ----------------
    est, bias, K0, K1 = run_kalman(acc, gyro, dt, args.Qa, args.Qb, args.R)
    comp = complementary(acc, gyro, dt, args.alpha)
    gyro_only = acc[0] + np.cumsum(gyro * dt)
    m = metrics(data, est, bias, args)

    print(f"參數: R={args.R:g}  Qa={args.Qa:g}  Qb={args.Qb:g}")
    print(f"K[0] 收斂值 {K0[-1]:.4f}   K[1] 收斂值 {K1[-1]:.4f}")
    print(f"靜置雜訊 std {m['static_std']:.4f} deg   bias 最大偏移 {m['bias_max']:.3f} deg/s")
    if "rms_true" in m:
        rms = lambda x: np.sqrt(np.mean((x - data["true"]) ** 2))
        print(f"RMS 誤差  加速度計 {rms(acc):.2f}  純陀螺儀 {rms(gyro_only):.2f}"
              f"  互補 {rms(comp):.2f}  卡爾曼 {m['rms_true']:.2f}  (deg)")
    print(f"動完後的靜止 {m['events']} 段,其中 {m['recovered']} 段在靜止期間恢復"
          f"  恢復時間 平均 {fmt(m['rec_mean'], ' s', 0, 2)} 最大 {fmt(m['rec_max'], ' s', 0, 2)}"
          f"  靜止殘差 RMS {fmt(m['resid_rms'], ' deg', 0)}")
    same_as_fw = (args.R, args.Qa, args.Qb) == (FW_R, FW_Q_ANGLE, FW_Q_BIAS)
    if data["mcu"] is not None and same_as_fw:
        diff = est - data["mcu"]
        print(f"Python 與 MCU 結果差異 RMS {np.sqrt(np.mean(diff**2)):.4f} deg (驗證 MCU 實作)")

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    has_true = data["true"] is not None
    n = 4
    fig, ax = plt.subplots(n, 1, figsize=(11, 3 * n), sharex=True)
    a0 = ax[0]
    a0.plot(t, acc, ".", ms=2, alpha=0.4, label="accel angle")
    a0.plot(t, gyro_only, label="gyro only")
    a0.plot(t, comp, label=f"complementary (a={args.alpha})")
    a0.plot(t, est, lw=2, label=f"Kalman R={args.R:g} Qb={args.Qb:g}")
    if data["mcu"] is not None and same_as_fw:
        a0.plot(t, data["mcu"], "--", lw=1, label="Kalman (MCU)")
    if has_true:
        a0.plot(t, data["true"], "k--", lw=1, label="true")
    lo = min(np.percentile(est, 1), np.percentile(gyro_only, 1)) - 10
    hi = max(np.percentile(est, 99), np.percentile(gyro_only, 99)) + 10
    a0.set_ylim(lo, hi)
    a0.set_ylabel("angle (deg)")
    a0.legend(loc="upper right", fontsize=8)

    still = m["_still"]
    for a_ in (ax[0], ax[1]):                      # 灰底 = 偵測為靜止的時段
        a_.fill_between(t, 0, 1, where=still, transform=a_.get_xaxis_transform(),
                        color="gray", alpha=0.12, lw=0)
    i = 1
    if has_true:
        ax[i].plot(t, est - data["true"], label="Kalman error (vs true)")
        ax[i].plot(t, comp - data["true"], alpha=0.7, label="complementary error")
    else:
        e_still = np.where(still, est - acc, np.nan)
        c_still = np.where(still, comp - acc, np.nan)
        ax[i].plot(t, e_still, label="Kalman - accel (still only)")
        ax[i].plot(t, c_still, alpha=0.7, label="complementary - accel (still only)")
        ax[i].axhline(args.tol, color="k", ls=":", lw=1)
        ax[i].axhline(-args.tol, color="k", ls=":", lw=1, label=f"±tol {args.tol} deg")
    ax[i].set_ylabel("error (deg)")
    ax[i].legend(fontsize=8)
    i += 1
    ax[i].plot(t, bias, label="estimated bias")
    if data["true_bias"] is not None:
        ax[i].plot(t, data["true_bias"], "k--", label="true bias")
    ax[i].set_ylabel("bias (deg/s)")
    ax[i].legend(fontsize=8)
    ax[i + 1].plot(t, K0, label="K[0]")
    ax[i + 1].plot(t, K1, label="K[1]")
    ax[i + 1].set_ylabel("gain")
    ax[i + 1].set_xlabel("time (s)")
    ax[i + 1].legend(fontsize=8)
    for a in ax:
        a.grid(True, alpha=0.3)

    stem = Path(args.csv).stem if args.csv else "sim"
    out = f"{stem}_R{args.R:g}_Qb{args.Qb:g}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print(f"已輸出 {out}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
