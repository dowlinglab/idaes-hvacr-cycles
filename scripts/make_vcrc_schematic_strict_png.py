import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch, Polygon
import numpy as np

plt.rcParams['font.family'] = 'DejaVu Sans'
fig, ax = plt.subplots(figsize=(14, 8), dpi=300)
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)
ax.axis('off')

# Colors
BLACK = 'black'
RED_FILL = '#f7d9dc'
BLUE_FILL = '#d9ecff'
YELLOW = '#ffe66d'
HEAT_RED = '#d62728'
HEAT_BLUE = '#1f77b4'
GRAY = '#7a7a7a'

# Strict rectangular loop
xL, xR = 2.2, 11.2
yB, yT = 1.7, 6.3
pipe_lw = 3.0

ax.plot([xL, xR], [yT, yT], color=BLACK, lw=pipe_lw)
ax.plot([xR, xR], [yT, yB], color=BLACK, lw=pipe_lw)
ax.plot([xR, xL], [yB, yB], color=BLACK, lw=pipe_lw)
ax.plot([xL, xL], [yB, yT], color=BLACK, lw=pipe_lw)

# Components (uniform spacing around loop)
cx = (xL + xR) / 2
cy = (yB + yT) / 2

# Condenser / Evaporator centered on top/bottom pipes and same width
hx_w, hx_h = 3.6, 1.0
cond_x, cond_y = cx - hx_w/2, yT - hx_h/2
evap_x, evap_y = cx - hx_w/2, yB - hx_h/2

cond = Rectangle((cond_x, cond_y), hx_w, hx_h, facecolor=RED_FILL, edgecolor=BLACK, lw=1.8)
evap = Rectangle((evap_x, evap_y), hx_w, hx_h, facecolor=BLUE_FILL, edgecolor=BLACK, lw=1.8)
ax.add_patch(cond); ax.add_patch(evap)

# Zig-zag coils centered inside HX (not touching edges)
def draw_coil(x, y, w, h, n=9):
    xs = np.linspace(x + 0.45, x + w - 0.45, n)
    ys = np.array([y + h*0.34 if i % 2 == 0 else y + h*0.66 for i in range(n)])
    ax.plot(xs, ys, color=BLACK, lw=1.4)

draw_coil(cond_x, cond_y, hx_w, hx_h)
draw_coil(evap_x, evap_y, hx_w, hx_h)

# Compressor (right side center, medium size)
comp_r = 0.52
comp_cx, comp_cy = xR, cy
ax.add_patch(Circle((comp_cx, comp_cy), comp_r, facecolor='white', edgecolor=BLACK, lw=1.8))

# Expansion valve inline on left side
vx, vy = xL, cy
ax.add_patch(Polygon([[vx-0.42, vy+0.28], [vx-0.42, vy-0.28], [vx, vy]], closed=True, facecolor='white', edgecolor=BLACK, lw=1.8))
ax.add_patch(Polygon([[vx+0.42, vy+0.28], [vx+0.42, vy-0.28], [vx, vy]], closed=True, facecolor='white', edgecolor=BLACK, lw=1.8))

# Uniform directional arrows ON pipes
arr_lw = 2.2
arr_ms = 16
arrows = [
    ((xR, yB+0.9), (xR, yT-0.9)),   # 1->2
    ((xR-1.4, yT), (xL+1.4, yT)),   # 2->3
    ((xL, yT-0.9), (xL, yB+0.9)),   # 3->4
    ((xL+1.4, yB), (xR-1.4, yB)),   # 4->1
]
for a, b in arrows:
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=arr_ms, lw=arr_lw, color=BLACK))

# Labels
ax.text(comp_cx+0.75, comp_cy, 'Compressor', fontsize=12, va='center')
ax.text(cx, yT+0.03, 'Condenser', fontsize=12, ha='center', va='center', weight='bold')
ax.text(vx, vy-0.75, 'Expansion valve', fontsize=11, ha='center')
ax.text(cx, yB-0.03, 'Evaporator', fontsize=12, ha='center', va='center', weight='bold')

# Compressor work arrow and label
ax.add_patch(FancyArrowPatch((xR+1.45, cy), (xR+0.58, cy), arrowstyle='-|>', mutation_scale=arr_ms, lw=2.0, color=GRAY))
ax.text(xR+1.5, cy+0.28, r'$W_{comp}$', fontsize=11, color=GRAY)

# State markers
states = {
    1: (xR-0.38, yB-0.34),
    2: (xR-0.38, yT+0.34),
    3: (xL+0.38, yT+0.34),
    4: (xL+0.38, yB-0.34),
}
for s, (sx, sy) in states.items():
    ax.add_patch(Circle((sx, sy), 0.16, facecolor=YELLOW, edgecolor=BLACK, lw=1.1))
    ax.text(sx, sy, str(s), fontsize=10, ha='center', va='center', weight='bold')

# Phase labels offset from pipes
ax.text(xR+0.2, yB-0.34, 'Low-P vapor', fontsize=9, va='center')
ax.text(xR+0.2, yT+0.34, 'High-P vapor', fontsize=9, va='center')
ax.text(xL-0.2, yT+0.34, 'High-P liquid', fontsize=9, ha='right', va='center')
ax.text(xL-0.2, yB-0.34, 'Low-P liquid + vapor', fontsize=9, ha='right', va='center')

# Ambient and cold-space indicators
ax.add_patch(FancyArrowPatch((cx+2.2, yT+0.15), (cx+2.2, yT+0.95), arrowstyle='-|>', mutation_scale=14, lw=2.0, color=HEAT_RED))
ax.text(cx+2.2, yT+1.05, r'Ambient, $T_0$', fontsize=10, ha='center')

ax.add_patch(FancyArrowPatch((cx+2.2, yB-0.95), (cx+2.2, yB-0.15), arrowstyle='-|>', mutation_scale=14, lw=2.0, color=HEAT_BLUE))
ax.text(cx+2.2, yB-1.15, r'Cold space, $T_{CL}$', fontsize=10, ha='center')

# Equations outside loop
ax.text(xR+0.55, cy-0.75, r'$W_{comp}=\dot{m}(h_2-h_1)$', fontsize=10)
ax.text(cx, yT+1.0, r'$Q_{cond}=\dot{m}(h_2-h_3)$', fontsize=10, ha='center')
ax.text(xL-1.55, cy-0.1, r'$h_3=h_4$', fontsize=10)
ax.text(cx, yB-1.0, r'$Q_{evap}=\dot{m}(h_1-h_4)$', fontsize=10, ha='center')
ax.text(xR+0.55, yB-1.0, r'$COP=\dfrac{Q_{evap}}{W_{comp}}$', fontsize=11)

# Title
ax.text(7, 7.55, 'Single-Stage Vapor Compression Cycle Model', fontsize=18, weight='bold', ha='center')

out = '/Users/snarasi2/idaes-hvacr-cycles/vcrc_single_stage_schematic_strict.png'
fig.savefig(out, bbox_inches='tight', dpi=300)
print(out)
