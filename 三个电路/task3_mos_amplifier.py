# -*- coding: utf-8 -*-
"""
任务③ NMOS 共源极放大器 (Common-Source Amplifier)
=================================================
题目固定参数 (题卡):
  VDD = 5 V, Rg1 = 60 kΩ, Rg2 = 40 kΩ, Rd = 2 kΩ, Cb1 足够大
  NMOS 模型: K = 0.8 mA/V², V_th = 1 V, λ = 0.02 /V
  输入 Vi = 10 mV (峰值), 1 kHz 正弦波

本脚本完成:
  1. 完整电路图 / 直流通路图 / 小信号等效模型   -> 3 张原理图
  2. 直流工作点 (.op):  V_GS, I_D, V_DS        -> 与手算对比
  3. 瞬态分析 (.tran):  输入/输出波形 (反相放大) -> 实测电压增益
  4. 交流分析 (.ac):    1 kHz 处复增益, 作旁证
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
# 一、参数
# =====================================================================
VDD = 5.0
RG1 = 60e3
RG2 = 40e3
RD = 2e3
CB1 = 100e-6          # 100 µF, 1kHz 时容抗仅 1.59 Ω, 满足 "足够大"
K = 0.8e-3            # A/V²  (题给导电因子)
VTH = 1.0             # V     (题给阈值电压)
LAM = 0.02            # 1/V   (题给沟道长度调制系数)
VI_PEAK = 10e-3       # 10 mV
F_IN = 1000.0

results = {}

# =====================================================================
# 二、手算静态工作点 (Q 点) —— 详细推导
#   直流通路: 电容 Cb1 对直流相当于开路(断开), 输入源不起作用;
#             栅极电流≈0 (MOS 管栅极绝缘, 不取电流),
#             故栅极电压完全由 Rg1、Rg2 分压决定:
#   V_G = VDD·Rg2/(Rg1+Rg2) = 5×40/100 = 2 V ; 源极接地 V_S = 0
#   V_GS = V_G − V_S = 2 V
#   过驱动电压 V_ov = V_GS − V_th = 2−1 = 1 V
#
#   饱和区电流公式(Shichman–Hodges, 含沟道长度调制 λ):
#     I_D = K·V_ov²·(1 + λ·V_DS),   V_DS = VDD − I_D·Rd
#   解出 I_D (本例 0.032 来自 λ·Rd = 0.02×2000):
#     I_D = K·V_ov²·(1+λ·VDD) / (1 + K·V_ov²·λ·Rd)
# =====================================================================
VG = VDD * RG2 / (RG1 + RG2)
VGS = VG                       # 源极接地
VOV = VGS - VTH
ID_th = K * VOV**2 * (1 + LAM * VDD) / (1 + K * VOV**2 * LAM * RD)
VDS_th = VDD - ID_th * RD
print(f'[手算] V_G = {VG:.3f} V, V_GS = {VGS:.3f} V, V_ov = {VOV:.3f} V')
print(f'[手算] I_D = {ID_th*1e3:.4f} mA')
print(f'[手算] V_DS = {VDS_th:.4f} V')

# ---- 饱和区判断: V_DS ≥ V_GS − V_th (=V_ov) 即工作在饱和区 ----
sat_margin = VDS_th - VOV
print(f'[判断] V_DS − V_ov = {sat_margin:.3f} V > 0  -> 饱和区 ✓')

# =====================================================================
# 三、手算小信号参数
#   gm (跨导 transconductance): 栅压对漏极电流的控制能力,
#       单位 S (西门子), 数值小常用 mS (毫西)
#       gm = 2K·V_ov = 1.6 mS
#   ro (管子输出电阻): 沟道长度调制带来的等效内阻
#       ro = 1/(λ·I_D) ≈ 58.6 kΩ
#   Av (电压增益 voltage gain):
#       忽略 ro : Av = −gm·Rd = −3.2  (负号=反相)
#       计入 ro : Av = −gm·(Rd∥ro) ≈ −3.09
# =====================================================================
gm_th = 2 * K * VOV
ro_th = 1 / (LAM * ID_th)
AV_ideal = -gm_th * RD
AV_precise = -gm_th * (RD * ro_th / (RD + ro_th))
# 严格按 Shichman-Hodges 小信号定义 (SPICE 数值微分口径):
#   gm_eff = 2K·Vov·(1+λV_DS),  g_ds = K·Vov²·λ,  Av = −gm_eff/(g_ds+1/Rd)
gm_eff = 2 * K * VOV * (1 + LAM * VDS_th)
gds = K * VOV**2 * LAM
ro_ds = 1 / gds
AV_exact = -gm_eff / (gds + 1 / RD)
print(f'[手算] gm = {gm_th*1e3:.3f} mS, ro = {ro_th/1e3:.2f} kΩ')
print(f'[手算] Av = {AV_ideal:.2f} (忽略ro), {AV_precise:.2f} (教材口径含ro)')
print(f'[手算] 严格模型: gm_eff = {gm_eff*1e3:.3f} mS, 1/gds = {ro_ds/1e3:.2f} kΩ, Av = {AV_exact:.3f}')

results.update(dict(VG=VG, VGS_theory=VGS, ID_theory=ID_th, VDS_theory=VDS_th,
                    VOV=VOV, sat_margin=sat_margin,
                    gm_theory=gm_th, ro_theory=ro_th,
                    gm_eff=gm_eff, gds=gds, ro_ds=ro_ds,
                    Av_ideal=AV_ideal, Av_precise=AV_precise, Av_exact=AV_exact))

# =====================================================================
# 四、画三张原理图
# =====================================================================
def vresistor(ax, x, y1, y2, label, color='#c0392b', lab_dx=0.35):
    """竖直电阻 (锯齿)"""
    ax.plot([x, x], [y1, (y1+y2)/2+0.55], 'k', lw=2)
    ys = np.linspace((y1+y2)/2+0.55, (y1+y2)/2-0.55, 12)
    xs = x + 0.17 * np.array([0, 1, -1]*3 + [0, 0, 0])
    ax.plot(xs, ys, 'k', lw=2)
    ax.plot([x, x], [(y1+y2)/2-0.55, y2], 'k', lw=2)
    ax.text(x+lab_dx, (y1+y2)/2, label, fontsize=10, color=color, va='center')

def ground(ax, x, y=1.2):
    ax.plot([x-0.32, x+0.32], [y, y], 'k', lw=2)
    ax.plot([x-0.21, x+0.21], [y-0.15, y-0.15], 'k', lw=2)
    ax.plot([x-0.1, x+0.1], [y-0.3, y-0.3], 'k', lw=2)

def mosfet(ax, x=5.0, yd=4.3, ys=2.5):
    """画 NMOS 管符号: 沟道竖线在 x, 栅线在 x-0.7, 源极箭头朝里"""
    yc = (yd + ys) / 2
    # 漏、源引线
    ax.plot([x, x], [yd, yc+0.55], 'k', lw=2)
    ax.plot([x, x], [yc-0.55, ys], 'k', lw=2)
    # 沟道 (增强型画成粗竖线) 与栅线 (有间隙)
    ax.plot([x, x], [yc+0.55, yc-0.55], 'k', lw=3)
    gx = x - 0.75
    ax.plot([gx, gx], [yc+0.6, yc-0.6], 'k', lw=2.5)
    # 栅极节点引线
    ax.plot([gx-1.3, gx], [yc, yc], 'k', lw=2)
    ax.plot(gx-1.3, yc, 'ko', ms=5)
    # 源极箭头 (NMOS: 箭头指向沟道, 即朝左上方)
    ax.annotate('', xy=(x-0.12, yc-0.28), xytext=(x-0.55, yc-0.72),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2))
    # 衬底 B 引线: 右侧箭头指向沟道, 再向下接源极
    bx = x + 0.75
    ax.annotate('', xy=(x+0.1, yc-0.1), xytext=(bx, yc-0.1),
                arrowprops=dict(arrowstyle='-|>', color='k', lw=2))
    ax.plot([bx, bx], [yc-0.1, ys], 'k', lw=2)
    ax.plot([x, bx], [ys, ys], 'k', lw=2)
    return gx-1.3, yc   # 栅极节点坐标

def draw_full(path):
    fig, ax = plt.subplots(figsize=(8.2, 7))
    ax.set_xlim(0, 8.5); ax.set_ylim(0.5, 6.8); ax.axis('off')
    # 顶部 VDD 母线
    ax.plot([3, 6.3], [6.3, 6.3], 'k', lw=2)
    ax.plot(6.3, 6.3, 'ko', mfc='white', ms=9)
    ax.text(6.55, 6.3, '$V_{DD}$ (5 V)', fontsize=13, va='center')
    # Rg1, Rd
    vresistor(ax, 3, 6.3, 3.4, '$R_{g1}$=60 kΩ')
    vresistor(ax, 5, 6.3, 4.3, '$R_d$=2 kΩ', lab_dx=0.3)
    # Rg2
    vresistor(ax, 3, 3.4, 1.2, '$R_{g2}$=40 kΩ')
    ground(ax, 3)
    # 栅极横线到 MOS
    ax.plot([3, 2.95], [3.4, 3.4], 'k', lw=2)
    # MOS 管 (栅节点在 x=2.95,y=3.4)
    gx, gy = mosfet(ax)
    ax.plot([5, 5], [2.5, 1.2], 'k', lw=2)
    ground(ax, 5)
    ax.plot(3, 3.4, 'ko', ms=5)
    # Cb1 + 输入源
    ax.plot([1.55, 2.95], [3.4, 3.4], 'k', lw=2)
    ax.plot([0.65, 1.45], [3.6, 3.6], 'k', lw=3)
    ax.plot([0.65, 1.45], [3.2, 3.2], 'k', lw=3)
    ax.text(1.05, 3.95, '$C_{b1}$\n(足够大)', ha='center', fontsize=10, color='#8e44ad')
    ax.add_patch(plt.Circle((0.55, 2.55), 0.42, fill=False, color='k', lw=2))
    ax.plot([0.55, 0.55], [2.97, 3.2], 'k', lw=2)
    ax.plot([0.55, 0.65], [3.2, 3.2], 'k', lw=2)
    ax.plot([0.55, 0.55], [2.13, 1.2], 'k', lw=2)
    ax.text(0.55, 2.62, '~', ha='center', fontsize=16)
    ax.text(0.02, 3.05, '$v_i$\n10mV\n1kHz', fontsize=10)
    ground(ax, 0.55)
    # 输出端子
    ax.plot(5, 4.3, 'ko', ms=5)
    ax.plot([5, 6.6], [4.3, 4.3], 'k', lw=2)
    ax.plot(6.6, 4.3, 'ko', mfc='white', ms=9)
    ax.text(6.8, 4.3, '$v_O$', fontsize=13, va='center')
    # 引脚字母
    ax.text(3.15, 3.55, 'g', fontsize=13, style='italic')
    ax.text(5.15, 4.45, 'd', fontsize=13, style='italic')
    ax.text(5.15, 2.55, 's', fontsize=13, style='italic')
    ax.text(5.95, 3.0, 'B', fontsize=12)
    ax.text(3.35, 2.6, 'T (NMOS)', fontsize=11, color='#2c3e50')
    ax.set_title('图3-1  NMOS 共源极放大器完整电路 (按题卡)', fontsize=13)
    plt.tight_layout(); fig.savefig(path, dpi=130); plt.close()

def draw_dc(path):
    """直流通路: Cb1 开路 -> 输入部分移除, 其余不变"""
    fig, ax = plt.subplots(figsize=(7.2, 6.6))
    ax.set_xlim(1, 8); ax.set_ylim(0.5, 6.8); ax.axis('off')
    ax.plot([3, 6.3], [6.3, 6.3], 'k', lw=2)
    ax.text(6.45, 6.3, '$V_{DD}$ (5 V)', fontsize=13, va='center')
    vresistor(ax, 3, 6.3, 3.4, '$R_{g1}$=60 kΩ')
    vresistor(ax, 5, 6.3, 4.3, '$R_d$=2 kΩ', lab_dx=0.3)
    vresistor(ax, 3, 3.4, 1.2, '$R_{g2}$=40 kΩ')
    ground(ax, 3)
    ax.plot([3, 2.95], [3.4, 3.4], 'k', lw=2)
    mosfet(ax)
    ax.plot([5, 5], [2.5, 1.2], 'k', lw=2)
    ground(ax, 5)
    ax.plot(3, 3.4, 'ko', ms=5)
    # Cb1 位置画断开的电容 + 叉号, 表示直流开路
    ax.plot([2.95, 2.4], [3.4, 3.4], 'k', lw=2)
    ax.plot([2.2, 2.2], [3.2, 3.6], 'r', lw=2)
    ax.plot([2.0, 2.0], [3.2, 3.6], 'r', lw=2)
    ax.plot([1.75, 2.45], [3.65, 3.15], 'r', lw=1.5)
    ax.plot([1.75, 2.45], [3.15, 3.65], 'r', lw=1.5)
    ax.text(1.2, 3.95, '$C_{b1}$ 对直流开路\n(输入信号被隔断)', fontsize=10, color='red')
    ax.text(3.15, 3.55, 'g', fontsize=13, style='italic')
    ax.text(5.15, 4.45, 'd', fontsize=13, style='italic')
    ax.text(5.15, 2.55, 's', fontsize=13, style='italic')
    ax.set_title('图3-2  直流通路图 (求静态工作点时, 电容视为开路)', fontsize=12)
    plt.tight_layout(); fig.savefig(path, dpi=130); plt.close()

def draw_small_signal(path):
    """小信号等效模型"""
    fig, ax = plt.subplots(figsize=(9, 5.2))
    ax.set_xlim(0, 11); ax.set_ylim(0, 6); ax.axis('off')
    # ---- 输入回路 ----
    ax.add_patch(plt.Circle((1, 3), 0.42, fill=False, color='k', lw=2))
    ax.text(1, 3.08, '+', ha='center', fontsize=12); ax.text(1, 2.55, '−', ha='center', fontsize=12)
    ax.text(1, 3.85, '$v_i$', ha='center', fontsize=13)
    ax.plot([1, 2.6], [3, 3], 'k', lw=2)
    ax.plot(2.6, 3, 'ko', ms=5)
    ax.text(2.7, 3.15, 'g', fontsize=13, style='italic')
    # Rg1//Rg2 等效偏置电阻到地
    ax.plot([2.6, 2.6], [3, 2.4], 'k', lw=2)
    ys = np.linspace(2.4, 1.5, 12)
    xs = 2.6 + 0.15*np.array([0,1,-1]*3+[0,0,0])
    ax.plot(xs, ys, 'k', lw=2)
    ax.plot([2.6, 2.6], [1.5, 1.1], 'k', lw=2)
    ground(ax, 2.6, y=1.1)
    ax.text(2.85, 1.95, '$R_{g1}∥R_{g2}$ = 24 kΩ', fontsize=10, color='#c0392b')
    # g-s 开路 (栅极绝缘)
    ax.plot([2.6, 2.6], [1.1, 0.75], 'w', lw=2)
    ax.plot([1, 1], [2.58, 1.1], 'k', lw=2)
    ax.plot([1, 2.6], [1.1, 1.1], 'k', lw=2)
    ax.text(1.7, 0.65, 's  (源极=交流地)', fontsize=10)
    ax.text(2.0, 2.2, '$v_{gs}$', fontsize=12, color='#1f618d')
    ax.annotate('', xy=(2.35, 2.9), xytext=(2.35, 1.25),
                arrowprops=dict(arrowstyle='<->', color='#1f618d', lw=1.2))
    # ---- 输出回路 ----
    # Rd 上端接交流地 (VDD 对交流短路到地)
    gx2 = 6.2
    ax.plot([gx2, gx2], [5.2, 4.5], 'k', lw=2)
    ys = np.linspace(4.5, 3.6, 12)
    xs = gx2 + 0.15*np.array([0,1,-1]*3+[0,0,0])
    ax.plot(xs, ys, 'k', lw=2)
    ax.plot([gx2, gx2], [3.6, 3.0], 'k', lw=2)
    ground(ax, gx2, y=5.2)
    ax.text(gx2+0.35, 5.2, '交流地 ($V_{DD}$ 对交流短路)', fontsize=9, va='center')
    ax.text(gx2-1.35, 4.05, '$R_d$=2 kΩ', fontsize=10, color='#c0392b')
    ax.plot(gx2, 3.0, 'ko', ms=5)
    # 压控电流源 gm vgs (圆圈+箭头) 从 d 流向地
    ax.plot([gx2, gx2], [3.0, 2.55], 'k', lw=2)
    ax.add_patch(plt.Circle((gx2, 2.1), 0.45, fill=False, color='k', lw=2))
    ax.annotate('', xy=(gx2, 1.75), xytext=(gx2, 2.45),
                arrowprops=dict(arrowstyle='-|>', color='#d35400', lw=2.2))
    ax.text(gx2-0.15, 2.1, '$g_m v_{gs}$', ha='right', fontsize=11, color='#d35400')
    ax.plot([gx2, gx2], [1.65, 1.1], 'k', lw=2)
    # ro 并联
    ax.plot([gx2+1.0, gx2+1.0], [3.0, 2.75], 'k', lw=2)
    ys = np.linspace(2.75, 1.55, 12)
    xs = gx2+1.0 + 0.15*np.array([0,1,-1]*3+[0,0,0])
    ax.plot(xs, ys, 'k', lw=2)
    ax.plot([gx2+1.0, gx2+1.0], [1.55, 1.1], 'k', lw=2)
    ground(ax, gx2+1.0, y=1.1)
    ax.plot([gx2, gx2+1.0], [3.0, 3.0], 'k', lw=2)
    ax.plot([gx2, gx2+1.0], [1.1, 1.1], 'k', lw=2)
    ax.text(gx2+1.25, 2.15, '$r_o$≈62.5 kΩ', fontsize=10, color='#c0392b')
    # vo 端子
    ax.plot([gx2, 9.2], [3.0, 3.0], 'k', lw=2)
    ax.plot(9.2, 3.0, 'ko', mfc='white', ms=9)
    ax.text(9.4, 3.0, '$v_o$', fontsize=13, va='center')
    ax.text(gx2+0.15, 3.15, 'd', fontsize=13, style='italic')
    ax.set_title('图3-3  低频小信号等效模型 ($v_o = -g_m v_{gs}(R_d∥r_o)$, 反相放大)', fontsize=12)
    plt.tight_layout(); fig.savefig(path, dpi=130); plt.close()

draw_full('circuit_3.png')
draw_dc('dc_path_3.png')
draw_small_signal('small_signal_3.png')

# =====================================================================
# 五、搭电路 (SPICE LEVEL=1 模型)
#   SPICE 公式: ID = 0.5·KP·(W/L)·Vov²·(1+λ·VDS)
#   令 W/L = 1, 则需 0.5·KP = K -> KP = 1.6 mA/V²
# =====================================================================
def build_circuit(ac=False):
    cct = Circuit('CS-amplifier')
    cct.model('nch', 'NMOS', level=1, Vto=VTH, kp=1.6e-3, lambda_=LAM,
              gamma=0.0, phi=0.6, is_=1e-15)
    cct.V('dd', 'vdd', cct.gnd, VDD)
    cct.R('g1', 'vdd', 'g', RG1)
    cct.R('g2', 'g', cct.gnd, RG2)
    cct.R('d', 'vdd', 'd', RD)
    # M(name, 漏d, 栅g, 源s, 衬底b)  —— 衬底与源极同接 GND
    cct.M(1, 'd', 'g', cct.gnd, cct.gnd, model='nch', w=1e-6, l=1e-6)
    if ac:
        cct.SinusoidalVoltageSource('i', 'vi', cct.gnd, amplitude=0,
                                    ac_magnitude=VI_PEAK)
    else:
        cct.SinusoidalVoltageSource('i', 'vi', cct.gnd, amplitude=VI_PEAK,
                                    offset=0, frequency=F_IN)
    cct.C(1, 'vi', 'g', CB1)
    return cct

# =====================================================================
# 六、直流工作点仿真
# =====================================================================
op = build_circuit().simulator(temperature=27, nominal_temperature=27).operating_point()
VG_sim = float(op.nodes['g'].as_ndarray()[0])
VS_sim = 0.0
VD_sim = float(op.nodes['d'].as_ndarray()[0])
VGS_sim = VG_sim - VS_sim
VDS_sim = VD_sim - VS_sim
ID_sim = (VDD - VD_sim) / RD          # 漏极电流 = Rd 上电流
print(f'[仿真] V_GS = {VGS_sim:.4f} V, I_D = {ID_sim*1e3:.4f} mA, V_DS = {VDS_sim:.4f} V')

results.update(dict(VGS_sim=VGS_sim, ID_sim=ID_sim, VDS_sim=VDS_sim,
                    VGS_err=abs(VGS_sim-VGS)/VGS*100,
                    ID_err=abs(ID_sim-ID_th)/ID_th*100,
                    VDS_err=abs(VDS_sim-VDS_th)/VDS_th*100))

# =====================================================================
# 七、瞬态仿真: 5 个周期, 实测反相电压增益
# =====================================================================
tr = build_circuit().simulator(temperature=27, nominal_temperature=27).transient(
    step_time=1e-6, end_time=5e-3)
t = np.array(tr.time.as_ndarray())
vi = np.array(tr.nodes['vi'].as_ndarray())
vg = np.array(tr.nodes['g'].as_ndarray())
vd = np.array(tr.nodes['d'].as_ndarray())

# 取最后 2 个周期(已达稳态), 去掉直流偏置后比较峰峰值
steady = t >= 3e-3
vi_ac = vi[steady] - np.mean(vi[steady])
vo_ac = vd[steady] - np.mean(vd[steady])
Av_tran = -(np.ptp(vo_ac) / np.ptp(vi_ac))     # 负号: 反相
print(f'[瞬态] 实测增益 |Av| = {abs(Av_tran):.3f} (反相), 输出摆幅 = {np.ptp(vo_ac)*1e3:.2f} mVpp')
results['Av_tran'] = Av_tran
results['vout_pp_mV'] = np.ptp(vo_ac)*1e3

fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 6.6), sharex=True)
a1.plot(t*1e3, vi*1e3, color='#1f77b4', lw=1.3, label='$v_i$ 输入 (10 mV 峰值)')
a1.set(ylabel='输入电压 (mV)', title='图3-4  输入信号')
a1.grid(alpha=0.3); a1.legend(loc='upper right')
a2.plot(t*1e3, vd, color='#d62728', lw=1.4, label='$v_D$ 漏极输出(含直流偏置)')
a2.axhline(VD_sim, color='gray', ls='--', lw=1, label=f'静态 $V_D$ = {VD_sim:.2f} V')
a2.set(xlabel='时间 t (ms)', ylabel='输出电压 (V)',
       title=f'图3-5  输出信号: 与输入反相, 峰峰值 {np.ptp(vo_ac)*1e3:.1f} mV (放大约 {abs(Av_tran):.1f} 倍)')
a2.grid(alpha=0.3); a2.legend(loc='upper right')
plt.tight_layout(); plt.savefig('waveforms_3.png', dpi=130); plt.close()

# 另画一张去偏置的输入/输出同框对比图(反相关系更直观)
fig, ax = plt.subplots(figsize=(9.5, 3.8))
ax.plot(t[steady]*1e3, vi_ac*1e3, color='#1f77b4', lw=1.4, label='$v_i$ (输入, mV)')
ax.plot(t[steady]*1e3, vo_ac*1e3, color='#d62728', lw=1.4,
        label=f'$v_o$ 交流分量 (mV, 幅度 ×{abs(Av_tran):.1f})')
ax.axhline(0, color='k', lw=0.8)
ax.set(xlabel='时间 t (ms)', ylabel='交流电压 (mV)',
       title='图3-6  去直流后输入 vs 输出: 相位相反 (反相放大器)')
ax.grid(alpha=0.3); ax.legend(loc='upper right')
plt.tight_layout(); plt.savefig('waveforms_ac_3.png', dpi=130); plt.close()

# =====================================================================
# 八、AC 分析旁证: 1 kHz 处复增益
# =====================================================================
ac = build_circuit(ac=True).simulator(temperature=27, nominal_temperature=27).ac(
    start_frequency=10, stop_frequency=1e6, number_of_points=120, variation='dec')
fa = np.array(ac.frequency.as_ndarray())
ha = np.array(ac.nodes['d'].as_ndarray())
i1k = np.argmin(np.abs(fa - F_IN))
Av_ac = ha[i1k] / VI_PEAK
print(f'[AC] 1kHz 增益 = {Av_ac.real:.3f} {Av_ac.imag:+.3f}j, |Av|={abs(Av_ac):.3f}, 相位={np.angle(Av_ac,deg=True):.1f}°')
results['Av_ac'] = complex(Av_ac.real, Av_ac.imag).__repr__()
results['Av_ac_mag'] = float(abs(Av_ac))
results['Av_ac_phase_deg'] = float(np.angle(Av_ac, deg=True))
results['gm_sim'] = abs(Av_ac) / (RD * (ro_th)/(RD+ro_th))  # 由增益反推 gm(近似)
results['Av_exact_err'] = abs(abs(Av_tran) - abs(AV_exact)) / abs(AV_exact) * 100

fig, (b1, b2) = plt.subplots(2, 1, figsize=(9, 5.8), sharex=True)
b1.semilogx(fa, 20*np.log10(np.abs(ha)/VI_PEAK), color='#c0392b', lw=1.6)
b1.axvline(F_IN, color='gray', ls='--', lw=1)
b1.set(ylabel='增益 (dB)', title='图3-7  共源放大器频响 (中频增益 ≈ %.1f dB ≈ %.2f 倍)'
        % (20*np.log10(abs(Av_ac)), abs(Av_ac)))
b1.grid(alpha=0.3, which='both')
b2.semilogx(fa, np.angle(ha, deg=True), color='#2980b9', lw=1.6)
b2.axvline(F_IN, color='gray', ls='--', lw=1)
b2.set(xlabel='频率 (Hz)', ylabel='相位 (°)')
b2.grid(alpha=0.3, which='both')
plt.tight_layout(); plt.savefig('bode_3.png', dpi=130); plt.close()

with open('results_3.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2, default=str)
print('\n任务③ 完成, 已生成 5 张图 + results_3.json')
