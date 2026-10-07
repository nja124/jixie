"""
RC 低通滤波器瞬态分析示例

电路:
    Vin ---[ R ]---+--- Vout
                   |
                  [ C ]
                   |
                  GND

输入: 方波 (PULSE)
分析: 瞬态分析 (.tran)
输出: 输入/输出波形对比图 + 理论 RC 时间常数
"""

import matplotlib
matplotlib.use('Agg')  # 非交互后端，保存为图片
import matplotlib.pyplot as plt
import numpy as np

# Windows 中文支持
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

from PySpice.Spice.Netlist import Circuit
import PySpice.Unit as U


def build_circuit(r_ohm, c_farad, v_low, v_high, period, rise_time=1e-9, fall_time=1e-9):
    """构建 RC 低通滤波器电路, 输入为方波"""
    circuit = Circuit('RC Low-Pass Filter')
    # PULSE 源: v1 v2 delay rise fall pulse_width period
    # 用 50% 占空比方波
    circuit.PulseVoltageSource('input', 'vin', circuit.gnd,
                               initial_value=v_low, pulsed_value=v_high,
                               pulse_width=period / 2, period=period,
                               rise_time=rise_time, fall_time=fall_time)
    circuit.R(1, 'vin', 'out', r_ohm)
    circuit.C(1, 'out', circuit.gnd, c_farad)
    return circuit


def main():
    # ----- 电路参数 -----
    R = 1.0 @ U.u_kΩ       # 1 kΩ
    C = 1.0 @ U.u_uF       # 1 µF
    V_LOW = 0.0            # 低电平 0 V
    V_HIGH = 5.0           # 高电平 5 V
    PERIOD = 10 @ U.u_ms   # 方波周期 10 ms (100 Hz)

    # RC 时间常数
    r_val = float(R)
    c_val = float(C)
    period_val = float(PERIOD)
    tau = r_val * c_val
    f_cutoff = 1.0 / (2.0 * np.pi * tau)
    print(f'R = {r_val:.0f} Ω')
    print(f'C = {c_val * 1e6:.1f} µF')
    print(f'τ = R·C = {tau * 1e3:.2f} ms')
    print(f'截止频率 fc = 1/(2πτ) = {f_cutoff:.1f} Hz')
    print(f'输入方波频率 = {1.0 / period_val:.0f} Hz')
    print()

    # ----- 构建并仿真 -----
    circuit = build_circuit(r_val, c_val, V_LOW, V_HIGH, float(PERIOD))
    simulator = circuit.simulator(temperature=25, nominal_temperature=25)

    # 瞬态分析: 步长 1µs, 仿真 5 个周期
    step = 1 @ U.u_us
    end_time = 5 * period_val
    analysis = simulator.transient(step_time=step, end_time=end_time)

    # ----- 提取波形 (numpy 2.x 兼容写法) -----
    time = np.array(analysis.time.as_ndarray())
    vin = np.array(analysis.nodes['vin'].as_ndarray())
    vout = np.array(analysis.nodes['out'].as_ndarray())

    # ----- 绘图 -----
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    # 上图: 完整波形
    ax1.plot(time * 1e3, vin, label='Vin (输入方波)', color='#1f77b4', linewidth=1.2)
    ax1.plot(time * 1e3, vout, label='Vout (RC 输出)', color='#d62728', linewidth=1.5)
    ax1.axhline(V_HIGH * (1 - np.exp(-1)), color='gray', linestyle='--', alpha=0.6,
                label=f'63.2% 充电电平 (τ={tau*1e3:.1f} ms)')
    ax1.set_ylabel('电压 (V)')
    ax1.set_title(f'RC 低通滤波器瞬态响应  (R={r_val:.0f}Ω, C={c_val*1e6:.1f}µF, τ={tau*1e3:.1f}ms, fc={f_cutoff:.0f}Hz)')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.3, V_HIGH + 0.5)

    # 下图: 放大第一个上升沿, 展示指数充电
    mask = (time >= 0) & (time <= period_val / 2)
    ax2.plot(time[mask] * 1e3, vin[mask], label='Vin', color='#1f77b4', linewidth=1.2)
    ax2.plot(time[mask] * 1e3, vout[mask], label='Vout', color='#d62728', linewidth=1.5)
    # 理论充电曲线: Vout(t) = V_HIGH * (1 - exp(-t/τ))
    t_theory = np.linspace(0, period_val / 2, 500)
    v_theory = V_HIGH * (1.0 - np.exp(-t_theory / tau))
    ax2.plot(t_theory * 1e3, v_theory, label='理论 1-e^(-t/τ)', color='green',
             linestyle=':', linewidth=2)
    ax2.set_xlabel('时间 (ms)')
    ax2.set_ylabel('电压 (V)')
    ax2.set_title('第一个上升沿: 指数充电过程')
    ax2.legend(loc='lower right')
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.2, V_HIGH + 0.2)

    plt.tight_layout()
    out_path = 'rc_lowpass_transient.png'
    plt.savefig(out_path, dpi=120)
    print(f'波形图已保存: {out_path}')

    # ----- 数值验证 -----
    # 在 t=τ 时, 仿真输出应接近 V_HIGH*(1-1/e) ≈ 0.632*V_HIGH
    idx_tau = np.argmin(np.abs(time - tau))
    vout_at_tau = vout[idx_tau]
    expected = V_HIGH * (1 - np.exp(-1))
    print(f'\n数值验证 (t=τ 时刻):')
    print(f'  仿真 Vout(τ) = {vout_at_tau:.4f} V')
    print(f'  理论 Vout(τ) = {expected:.4f} V  (63.2% × {V_HIGH}V)')
    print(f'  误差 = {abs(vout_at_tau - expected):.4f} V ({abs(vout_at_tau-expected)/expected*100:.2f}%)')


if __name__ == '__main__':
    main()
