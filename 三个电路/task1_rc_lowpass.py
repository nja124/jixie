# -*- coding: utf-8 -*-
"""
任务① RC 低通滤波器 (RC Low-Pass Filter)
=========================================
电路连接:
    Vin (方波/交流源) ----[ R ]----+---- Vout (输出)
                                   |
                                  [ C ]
                                   |
                                  GND (地, 0V)

本脚本完成 3 件事:
  1. 画出电路连接示意图        -> circuit_1.png
  2. 方波瞬态分析 (.tran)      -> transient_1.png
  3. 交流波特图 (.ac)          -> bode_1.png
同时把 "手算理论值" 和 "仿真测量值" 写入 results_1.json, 供报告使用。
"""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import u_kΩ, u_uF, u_ms, u_us, u_V, u_Hz

# =====================================================================
# 一、元件参数 (本题允许自定)
# =====================================================================
R = 1.0 @ u_kΩ          # 电阻 R = 1 千欧 = 1000 Ω
C = 1.0 @ u_uF          # 电容 C = 1 微法 = 1e-6 F
V_HIGH = 5.0            # 方波高电平 5 V
V_LOW = 0.0             # 方波低电平 0 V
F_SQUARE = 100.0        # 方波频率 100 Hz, 周期 T = 10 ms

r, c = float(R), float(C)

# =====================================================================
# 二、手算理论值
#   τ (tau, 时间常数): 电容充放电的快慢, τ = R·C, 单位秒(s)
#   fc (截止频率):     信号幅度衰减到 70.7% (功率一半, -3dB) 时的频率
#                      fc = 1 / (2πRC)
# =====================================================================
tau_th = r * c                       # 理论时间常数 (s)
fc_th = 1.0 / (2 * np.pi * tau_th)   # 理论截止频率 (Hz)
results = {
    'R_ohm': r, 'C_F': c,
    'tau_theory_ms': tau_th * 1e3,
    'fc_theory_hz': fc_th,
}
print(f'[手算] τ = R·C = {r:.0f} × {c:.2e} = {tau_th*1e3:.3f} ms')
print(f'[手算] fc = 1/(2πRC) = {fc_th:.2f} Hz')

# =====================================================================
# 三、画电路示意图
# =====================================================================
def draw_circuit(path):
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.axis('off')
    ax.set_xlim(0, 10); ax.set_ylim(0, 5)

    # 主导线: Vin -- R -- out
    ax.plot([1, 3], [3, 3], 'k', lw=2)                      # 源到R
    # 电阻 R (锯齿形)
    xs = np.linspace(3, 5, 14)
    ys = 3 + 0.22 * np.array([0, 1, -1] * 4 + [0, 0])
    ax.plot(xs, ys, 'k', lw=2)
    ax.plot([5, 7.2], [3, 3], 'k', lw=2)                    # R 到 out 节点
    ax.plot(7.2, 3, 'ko', ms=5)                             # out 节点
    ax.plot([7.2, 9], [3, 3], 'k', lw=2)                    # 输出端子引线
    ax.plot(9, 3, 'ko', mfc='white', ms=9)                  # 输出端子(空心)
    ax.text(9.25, 3, '$V_{out}$', va='center', fontsize=13)

    # 输入源 (圆圈内画正负)
    ax.add_patch(plt.Circle((1, 3), 0.45, fill=False, color='k', lw=2))
    ax.text(1, 3.12, '+', ha='center', fontsize=13)
    ax.text(1, 2.55, '−', ha='center', fontsize=13)
    ax.text(1, 3.85, '$V_{in}$', ha='center', fontsize=13)
    # 源接地
    ax.plot([1, 1], [2.55, 1.3], 'k', lw=2)

    # 电容 C (两块平行板) 从 out 节点向下
    ax.plot([7.2, 7.2], [3, 2.25], 'k', lw=2)
    ax.plot([6.85, 7.55], [2.25, 2.25], 'k', lw=3)
    ax.plot([6.85, 7.55], [1.95, 1.95], 'k', lw=3)
    ax.plot([7.2, 7.2], [1.95, 1.3], 'k', lw=2)

    # 接地符号 (三条渐短横线)
    ax.plot([1, 1], [1.3, 1.3], 'k', lw=2)
    for x0 in (1, 7.2):
        ax.plot([x0-0.4, x0+0.4], [1.3, 1.3], 'k', lw=2)
        ax.plot([x0-0.27, x0+0.27], [1.12, 1.12], 'k', lw=2)
        ax.plot([x0-0.14, x0+0.14], [0.94, 0.94], 'k', lw=2)

    ax.text(4, 3.75, 'R = 1 kΩ', ha='center', fontsize=12, color='#c0392b')
    ax.text(7.95, 2.1, 'C = 1 µF', fontsize=12, color='#c0392b')
    ax.text(0.15, 3, '输入', fontsize=11, color='#2c3e50')
    ax.set_title('图1-1  RC 低通滤波器电路示意图', fontsize=13)
    plt.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)

draw_circuit('circuit_1.png')

# =====================================================================
# 四、瞬态仿真: 输入 100Hz 方波, 观察电容充放电
# =====================================================================
ct = Circuit('RC-transient')
ct.PulseVoltageSource('in', 'vin', ct.gnd,
                      initial_value=V_LOW, pulsed_value=V_HIGH,
                      delay_time=0, rise_time=1e-9, fall_time=1e-9,
                      pulse_width=10e-3 / 2, period=10e-3)
ct.R(1, 'vin', 'out', r)
ct.C(1, 'out', ct.gnd, c)

sim = ct.simulator(temperature=25, nominal_temperature=25)
tr = sim.transient(step_time=1 @ u_us, end_time=50 @ u_ms)
t = np.array(tr.time.as_ndarray())
vin_t = np.array(tr.nodes['vin'].as_ndarray())
vout_t = np.array(tr.nodes['out'].as_ndarray())

# ---- 从充电曲线拟合仿真 τ: Vout(t) = V_H (1 - exp(-t/τ)) ----
half = (t >= 0) & (t <= 5e-3)
t_fit, v_fit = t[half], vout_t[half]
def charge(t, tau, vmax):
    return vmax * (1 - np.exp(-t / tau))
popt, _ = curve_fit(charge, t_fit, v_fit, p0=[1e-3, 5.0])
tau_sim = float(popt[0])
print(f'[仿真] 拟合 τ = {tau_sim*1e3:.3f} ms')
results['tau_sim_ms'] = tau_sim * 1e3

# ---- 瞬态波形图 (全貌 + 第一个上升沿放大) ----
fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6.5))
a1.plot(t*1e3, vin_t, color='#1f77b4', lw=1.1, label='Vin 输入方波')
a1.plot(t*1e3, vout_t, color='#d62728', lw=1.6, label='Vout 电容两端电压')
a1.axhline(V_HIGH*(1-np.exp(-1)), color='gray', ls='--', lw=1,
           label=f'63.2% 电平 = {V_HIGH*(1-np.exp(-1)):.2f} V (t = τ)')
a1.set(xlabel='时间 t (ms)', ylabel='电压 (V)', ylim=(-0.3, 5.6),
       title='图1-2  方波输入下的瞬态波形 (100 Hz, 5 个周期)')
a1.legend(fontsize=9); a1.grid(alpha=0.3)

a2.plot(t_fit*1e3, v_fit, color='#d62728', lw=1.8, label='仿真 Vout')
tt = np.linspace(0, 5e-3, 300)
a2.plot(tt*1e3, V_HIGH*(1-np.exp(-tt/tau_sim)), '--', color='green', lw=1.5,
        label='理论曲线 5(1−e$^{-t/τ}$)')
a2.axvline(tau_sim*1e3, color='purple', ls=':', lw=1.2)
a2.plot(tau_sim*1e3, V_HIGH*(1-np.exp(-1)), 'o', color='purple')
a2.annotate(f't = τ = {tau_sim*1e3:.2f} ms\nV = 3.16 V (63.2%)',
            xy=(tau_sim*1e3, V_HIGH*(1-np.exp(-1))), xytext=(2.3, 1.6),
            arrowprops=dict(arrowstyle='->', color='purple'), fontsize=10, color='purple')
a2.set(xlabel='时间 t (ms)', ylabel='电压 (V)', ylim=(0, 5.4),
       title='图1-3  第一个上升沿: 指数充电过程')
a2.legend(fontsize=9); a2.grid(alpha=0.3)
plt.tight_layout(); plt.savefig('transient_1.png', dpi=130); plt.close()

# =====================================================================
# 五、交流分析 (.ac): 扫频得到波特图, 找 -3dB 截止频率
# =====================================================================
ca = Circuit('RC-ac')
# ac_magnitude=1 表示交流小信号幅度 1V (仅用于 .ac 分析)
ca.SinusoidalVoltageSource('in', 'vin', ca.gnd, amplitude=0, ac_magnitude=1)
ca.R(1, 'vin', 'out', r)
ca.C(1, 'out', ca.gnd, c)
sima = ca.simulator(temperature=25, nominal_temperature=25)
ac = sima.ac(start_frequency=1.0, stop_frequency=100e3,
             number_of_points=200, variation='dec')
freq = np.array(ac.frequency.as_ndarray())
h = np.array(ac.nodes['out'].as_ndarray())          # Vout (Vin=1)
mag_db = 20 * np.log10(np.abs(h))                  # 幅度, 单位 dB
phase = np.angle(h, deg=True)                      # 相位, 单位度

# -3dB 截止频率: |H| = 1/√2 即 -3.01 dB
target_db = 20*np.log10(1/np.sqrt(2))
# 在单调下降段线性插值
sel = (freq > 10) & (freq < 5000)
fc_sim = float(np.interp(target_db, mag_db[sel][::-1], np.log10(freq[sel])[::-1]))
fc_sim = 10 ** fc_sim
print(f'[仿真] -3dB 截止频率 fc = {fc_sim:.2f} Hz')
results['fc_sim_hz'] = fc_sim

fig, (b1, b2) = plt.subplots(2, 1, figsize=(9, 6.2), sharex=True)
b1.semilogx(freq, mag_db, color='#c0392b', lw=1.8)
b1.axhline(-3.01, color='gray', ls='--', lw=1)
b1.axvline(fc_th, color='green', ls=':', lw=1.5, label=f'手算 fc = {fc_th:.1f} Hz')
b1.axvline(fc_sim, color='purple', ls=':', lw=1.5, label=f'仿真 fc = {fc_sim:.1f} Hz')
b1.set(ylabel='增益 20lg|H| (dB)', title='图1-4  波特图 (Bode Plot) — 幅频特性')
b1.legend(fontsize=9); b1.grid(alpha=0.3, which='both'); b1.set_ylim(-45, 3)
b2.semilogx(freq, phase, color='#2980b9', lw=1.8)
b2.axvline(fc_th, color='green', ls=':', lw=1.5)
b2.axhline(-45, color='gray', ls='--', lw=1)
b2.set(xlabel='频率 f (Hz, 对数坐标)', ylabel='相位 φ (°)',
       title='相频特性 (fc 处相位滞后 45°)')
b2.grid(alpha=0.3, which='both'); b2.set_ylim(-95, 5)
plt.tight_layout(); plt.savefig('bode_1.png', dpi=130); plt.close()

# 误差百分比
results['tau_err_pct'] = abs(tau_sim - tau_th) / tau_th * 100
results['fc_err_pct'] = abs(fc_sim - fc_th) / fc_th * 100
with open('results_1.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print('\n任务① 完成, 已生成: circuit_1.png, transient_1.png, bode_1.png, results_1.json')
