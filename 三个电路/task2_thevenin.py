﻿# -*- coding: utf-8 -*-
"""
任务② 戴维南定理 (Thevenin's Theorem) 验证
===========================================
定理内容 (大白话):
  任何一个 "含电源 + 一堆电阻" 的二端网络, 从两个接线端子 (A、B) 看进去,
  都可以等效替换成:  一个理想电压源 V_th 串联一个电阻 R_th。
     - V_th = 端口开路电压 V_oc   (什么都不接时 A、B 之间的电压)
     - R_th = 电源置零后端口等效电阻 = V_oc / I_sc
       (I_sc = 把 A、B 直接短接时流过短接线的电流)

本题自定参数:
  V1 = 12 V, R1 = 4 kΩ, R2 = 8 kΩ, 外接负载 RL = 4 kΩ

        V1(+12V)──[R1 4k]──┬── A (输出端口)
                            |
                           [R2 8k]
                            |
        GND ────────────────┴── B (B 与地同电位, 0V)

要跑的 3 次仿真:
  (a) 端口开路 -> 测 V_oc
  (b) 端口短路 -> 测 I_sc
  (c) 原电路 / 戴维南等效电路 分别接上同一个负载 RL -> 比较 V_L、I_L
"""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

from PySpice.Spice.Netlist import Circuit

# =====================================================================
# 一、元件参数
# =====================================================================
V1 = 12.0      # 直流电压源 12 V
R1 = 4e3       # 4 kΩ
R2 = 8e3       # 8 kΩ
RL = 4e3       # 负载电阻 4 kΩ

# =====================================================================
# 二、手算理论值
#   V_oc : R1、R2 串联分压  V_oc = V1·R2/(R1+R2)
#   R_th : 电压源短路(用导线代替)后, R1、R2 变成并联  R_th = R1//R2
#   I_sc : 短接线电流 = V_oc/R_th  (此时 R2 被短接, 也等于 V1/R1)
#   接负载: V_L = V_oc·RL/(R_th+RL), I_L = V_L/RL
# =====================================================================
Voc_th = V1 * R2 / (R1 + R2)
Rth_th = R1 * R2 / (R1 + R2)
Isc_th = Voc_th / Rth_th
VL_th = Voc_th * RL / (Rth_th + RL)
IL_th = VL_th / RL
print(f'[手算] V_oc = {Voc_th:.4f} V')
print(f'[手算] R_th = {Rth_th/1e3:.4f} kΩ')
print(f'[手算] I_sc = {Isc_th*1e3:.4f} mA')
print(f'[手算] 接负载: V_L = {VL_th:.4f} V, I_L = {IL_th*1e3:.4f} mA')

results = {
    'V1': V1, 'R1': R1, 'R2': R2, 'RL': RL,
    'Voc_theory': Voc_th, 'Rth_theory': Rth_th, 'Isc_theory': Isc_th,
    'VL_theory': VL_th, 'IL_theory': IL_th,
}

def run_op(circuit):
    """执行一次直流工作点仿真 (.op)"""
    sim = circuit.simulator(temperature=25, nominal_temperature=25)
    return sim.operating_point()

# =====================================================================
# 三、仿真 (a): 端口开路 —— 什么都不接, 测 A 点电压
# =====================================================================
c_oc = Circuit('open-circuit')
c_oc.V(1, 'vs', c_oc.gnd, V1)          # 12V 源: +接 vs, -接地
c_oc.R(1, 'vs', 'A', R1)
c_oc.R(2, 'A', c_oc.gnd, R2)
op = run_op(c_oc)
Voc_sim = float(op.nodes['A'].as_ndarray()[0])
print(f'[仿真] V_oc = {Voc_sim:.4f} V')
results['Voc_sim'] = Voc_sim

# =====================================================================
# 四、仿真 (b): 端口短路 —— A、B 之间接一根导线(0V 电压源当电流表),
#                测流过它的电流
# =====================================================================
c_sc = Circuit('short-circuit')
c_sc.V(1, 'vs', c_sc.gnd, V1)
c_sc.R(1, 'vs', 'A', R1)
c_sc.R(2, 'A', c_sc.gnd, R2)
c_sc.V('sc', 'A', c_sc.gnd, 0)          # 0V 源 = 导线, 同时充当电流表
op = run_op(c_sc)
# 流过电压源 V_sc 的电流 (SPICE 中支路电流按源名访问)
Isc_sim = None
for key in ('vsc', 'Vsc', 'V_sc'):
    try:
        Isc_sim = abs(float(op.branches[key].as_ndarray()[0]))
        break
    except Exception:
        continue
if Isc_sim is None:
    # 备用方法: 用 1 mΩ 小电阻代替短接线
    c_sc2 = Circuit('short-circuit-2')
    c_sc2.V(1, 'vs', c_sc2.gnd, V1)
    c_sc2.R(1, 'vs', 'A', R1)
    c_sc2.R(2, 'A', c_sc2.gnd, R2)
    c_sc2.R('sc', 'A', c_sc2.gnd, 1e-3)
    op2 = run_op(c_sc2)
    va = float(op2.nodes['A'].as_ndarray()[0])
    Isc_sim = va / 1e-3
print(f'[仿真] I_sc = {Isc_sim*1e3:.4f} mA')
results['Isc_sim'] = Isc_sim
Rth_sim = Voc_sim / Isc_sim
results['Rth_sim'] = Rth_sim
print(f'[仿真] R_th = V_oc/I_sc = {Rth_sim/1e3:.4f} kΩ')

# =====================================================================
# 五、仿真 (c1): 原电路接上负载 RL, 测 V_L、I_L
# =====================================================================
c_load = Circuit('original-with-load')
c_load.V(1, 'vs', c_load.gnd, V1)
c_load.R(1, 'vs', 'A', R1)
c_load.R(2, 'A', c_load.gnd, R2)
c_load.R('L', 'A', c_load.gnd, RL)
op = run_op(c_load)
VL_orig = float(op.nodes['A'].as_ndarray()[0])
IL_orig = VL_orig / RL
print(f'[仿真-原电路] V_L = {VL_orig:.4f} V, I_L = {IL_orig*1e3:.4f} mA')
results['VL_original'] = VL_orig
results['IL_original'] = IL_orig

# =====================================================================
# 六、仿真 (c2): 戴维南等效电路接上同一个负载
#   V_th (用仿真实测的 V_oc) 串联 R_th (用仿真实测值) 再接 RL
# =====================================================================
c_th = Circuit('thevenin-with-load')
c_th.V('th', 'vth', c_th.gnd, Voc_sim)   # 等效电压源
c_th.R('th', 'vth', 'A', Rth_sim)        # 等效内阻
c_th.R('L', 'A', c_th.gnd, RL)           # 同一个负载
op = run_op(c_th)
VL_equiv = float(op.nodes['A'].as_ndarray()[0])
IL_equiv = VL_equiv / RL
print(f'[仿真-等效电路] V_L = {VL_equiv:.4f} V, I_L = {IL_equiv*1e3:.4f} mA')
results['VL_equivalent'] = VL_equiv
results['IL_equivalent'] = IL_equiv

# =====================================================================
# 七、画电路图: 左=原网络(含负载), 右=戴维南等效电路
# =====================================================================
def zigzag(ax, x1, x2, y):
    xs = np.linspace(x1, x2, 12)
    ys = y + 0.18 * np.array([0, 1, -1] * 3 + [0, 0, 0])
    ax.plot(xs, ys, 'k', lw=2)

def ground(ax, x, y=1.0):
    ax.plot([x-0.35, x+0.35], [y, y], 'k', lw=2)
    ax.plot([x-0.23, x+0.23], [y-0.16, y-0.16], 'k', lw=2)
    ax.plot([x-0.11, x+0.11], [y-0.32, y-0.32], 'k', lw=2)

def draw_diagrams(path):
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.2))
    # ---- 左图: 原电路 ----
    axL.set_xlim(0, 8); axL.set_ylim(0, 6); axL.axis('off')
    axL.add_patch(plt.Circle((1, 4.5), 0.42, fill=False, color='k', lw=2))
    axL.text(1, 4.62, '+', ha='center', fontsize=12); axL.text(1, 4.1, '−', ha='center', fontsize=12)
    axL.text(1, 5.25, 'V1 = 12 V', ha='center', fontsize=11)
    axL.plot([1, 1], [4.08, 2.8], 'k', lw=2)
    axL.plot([1, 2.3], [2.8, 2.8], 'k', lw=2)
    zigzag(axL, 2.3, 3.7, 2.8)                       # R1
    axL.text(3, 3.35, 'R1 = 4 kΩ', ha='center', fontsize=10, color='#c0392b')
    axL.plot([3.7, 6.2], [2.8, 2.8], 'k', lw=2)
    axL.plot(6.2, 2.8, 'ko', ms=5)
    axL.plot(6.2, 2.8, 'ko', mfc='white', ms=9)
    axL.text(6.4, 2.95, 'A', fontsize=13, color='#1f618d', weight='bold')
    # R2 向下 (竖直电阻)
    axL.plot([5, 5], [2.8, 2.35], 'k', lw=2)
    ys = np.linspace(2.35, 1.15, 12)
    xs = 5 + 0.18 * np.array([0, 1, -1] * 3 + [0, 0, 0])
    axL.plot(xs, ys, 'k', lw=2)
    axL.plot([5, 5], [1.15, 1.0], 'k', lw=2)
    axL.text(5.35, 1.75, 'R2 = 8 kΩ', fontsize=10, color='#c0392b')
    ground(axL, 5)
    axL.plot([1, 1], [2.8, 1.0], 'k', lw=2)        # 源到地
    ground(axL, 1)
    # 负载 RL (虚线框表示外接)
    axL.plot([6.2, 7.3], [2.8, 2.8], 'k', lw=2)
    ys2 = np.linspace(2.8, 1.55, 12)
    xs2 = 7.3 + 0.18 * np.array([0, 1, -1] * 3 + [0, 0, 0])
    axL.plot(xs2, ys2, 'k', lw=2)
    axL.plot([7.3, 7.3], [1.55, 1.0], 'k', lw=2)
    ground(axL, 7.3)
    axL.text(7.55, 2.0, 'RL = 4 kΩ\n(负载)', fontsize=10, color='#27ae60')
    axL.plot(6.2, 1.0, 'ko', ms=4)
    axL.plot([6.2, 7.3], [1.0, 1.0], 'k', lw=1)
    axL.text(6.4, 0.6, 'B', fontsize=13, color='#1f618d', weight='bold')
    axL.set_title('图2-1  原含源二端网络 + 负载 (A、B 为端口)', fontsize=12)

    # ---- 右图: 戴维南等效 ----
    axR.set_xlim(0, 8); axR.set_ylim(0, 6); axR.axis('off')
    axR.add_patch(plt.Circle((1.5, 4.5), 0.42, fill=False, color='k', lw=2))
    axR.text(1.5, 4.62, '+', ha='center', fontsize=12); axR.text(1.5, 4.1, '−', ha='center', fontsize=12)
    axR.text(1.5, 5.25, 'V_th = V_oc = 8 V', ha='center', fontsize=11)
    axR.plot([1.5, 1.5], [4.08, 2.8], 'k', lw=2)
    zigzag(axR, 1.5, 3.6, 2.8)
    axR.text(2.55, 3.35, 'R_th = 2.67 kΩ', ha='center', fontsize=10, color='#c0392b')
    axR.plot([3.6, 6.2], [2.8, 2.8], 'k', lw=2)
    axR.plot(6.2, 2.8, 'ko', mfc='white', ms=9)
    axR.text(6.4, 2.95, 'A', fontsize=13, color='#1f618d', weight='bold')
    axR.plot([1.5, 1.5], [2.8, 1.0], 'k', lw=2)
    ground(axR, 1.5)
    ys2 = np.linspace(2.8, 1.55, 12)
    xs2 = 7.3 + 0.18 * np.array([0, 1, -1] * 3 + [0, 0, 0])
    axR.plot([6.2, 7.3], [2.8, 2.8], 'k', lw=2)
    axR.plot(xs2, ys2, 'k', lw=2)
    axR.plot([7.3, 7.3], [1.55, 1.0], 'k', lw=2)
    ground(axR, 7.3)
    axR.text(7.55, 2.0, 'RL = 4 kΩ\n(同一负载)', fontsize=10, color='#27ae60')
    axR.plot([6.2, 7.3], [1.0, 1.0], 'k', lw=1)
    axR.text(6.4, 0.6, 'B', fontsize=13, color='#1f618d', weight='bold')
    axR.set_title('图2-2  戴维南等效电路 + 同一负载', fontsize=12)

    plt.tight_layout(); fig.savefig(path, dpi=130); plt.close()

draw_diagrams('circuit_2.png')

# =====================================================================
# 八、结果汇总 (相对误差)
# =====================================================================
results['Voc_err'] = abs(Voc_sim - Voc_th) / Voc_th * 100
results['Rth_err'] = abs(Rth_sim - Rth_th) / Rth_th * 100
results['Isc_err'] = abs(Isc_sim - Isc_th) / Isc_th * 100
results['VL_diff'] = abs(VL_orig - VL_equiv)
results['IL_diff'] = abs(IL_orig - IL_equiv)
with open('results_2.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print('\n任务② 完成, 已生成: circuit_2.png, results_2.json')
