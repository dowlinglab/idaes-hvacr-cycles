import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch, Polygon
import numpy as np

plt.rcParams['font.family'] = 'DejaVu Sans'

fig, ax = plt.subplots(figsize=(14, 8), dpi=200)
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)
ax.axis('off')

# Colors
pipe_c = 'black'
cond_fill = '#f8d7da'  # light red/pink
evap_fill = '#d8ecff'  # light blue
state_fill = '#ffe66d' # yellow
heat_red = '#d62728'
heat_blue = '#1f77b4'
gray = '#7f7f7f'

# Main components positions
cond_x, cond_y, cond_w, cond_h = 4.0, 5.8, 6.0, 1.1
evap_x, evap_y, evap_w, evap_h = 4.0, 1.1, 6.0, 1.1
comp_cx, comp_cy, comp_r = 11.4, 3.95, 0.8

# Condenser and evaporator blocks
cond = Rectangle((cond_x, cond_y), cond_w, cond_h, facecolor=cond_fill, edgecolor='black', linewidth=2)
evap = Rectangle((evap_x, evap_y), evap_w, evap_h, facecolor=evap_fill, edgecolor='black', linewidth=2)
ax.add_patch(cond)
ax.add_patch(evap)

# Coil zig-zag inside heat exchangers
def add_coil(x0, y0, w, h, turns=11):
    xs = np.linspace(x0 + 0.35, x0 + w - 0.35, turns)
    ys = np.array([y0 + h*0.28 if i % 2 == 0 else y0 + h*0.72 for i in range(turns)])
    ax.plot(xs, ys, color='black', linewidth=1.8)

add_coil(cond_x, cond_y, cond_w, cond_h)
add_coil(evap_x, evap_y, evap_w, evap_h)

# Compressor symbol
comp = Circle((comp_cx, comp_cy), comp_r, facecolor='white', edgecolor='black', linewidth=2)
ax.add_patch(comp)
# compressor impeller hint
for ang in [30, 150, 270]:
    t = np.deg2rad(ang)
    ax.plot([comp_cx, comp_cx + 0.45*np.cos(t)], [comp_cy, comp_cy + 0.45*np.sin(t)], color='black', linewidth=1.4)

# Expansion valve symbol on left side (inline throttle)
valve_x = 2.6
valve_y = 3.95
tri1 = Polygon([[valve_x-0.45, valve_y+0.35],[valve_x-0.45, valve_y-0.35],[valve_x, valve_y]], closed=True, facecolor='white', edgecolor='black', linewidth=2)
tri2 = Polygon([[valve_x+0.45, valve_y+0.35],[valve_x+0.45, valve_y-0.35],[valve_x, valve_y]], closed=True, facecolor='white', edgecolor='black', linewidth=2)
ax.add_patch(tri1)
ax.add_patch(tri2)

# Pipe loop coordinates and lines (continuous rectangular loop with side devices)
# Top (2->3): condenser outlet/inlet points
xL = 2.6
xR = comp_cx
yT = 6.35
yB = 1.65

# segments
ax.plot([xR, cond_x+cond_w], [yT, yT], color=pipe_c, linewidth=2.5)        # short from comp to condenser right edge
ax.plot([cond_x, xL], [yT, yT], color=pipe_c, linewidth=2.5)                # condenser left to valve side
ax.plot([xL, xL], [yT, yB], color=pipe_c, linewidth=2.5)                    # left vertical through valve
ax.plot([xL, evap_x], [yB, yB], color=pipe_c, linewidth=2.5)                # valve to evaporator left
ax.plot([evap_x+evap_w, xR], [yB, yB], color=pipe_c, linewidth=2.5)         # evaporator right to compressor
ax.plot([xR, xR], [yB, yT], color=pipe_c, linewidth=2.5)                    # right vertical through compressor

# connect through compressor and exchangers visually
ax.plot([cond_x+cond_w, cond_x+cond_w], [yT, cond_y+cond_h/2], color=pipe_c, linewidth=2.5)
ax.plot([cond_x+cond_w, xR], [cond_y+cond_h/2, cond_y+cond_h/2], color=pipe_c, linewidth=2.5)
ax.plot([cond_x, cond_x], [yT, cond_y+cond_h/2], color=pipe_c, linewidth=2.5)
ax.plot([xL, cond_x], [cond_y+cond_h/2, cond_y+cond_h/2], color=pipe_c, linewidth=2.5)

ax.plot([evap_x, evap_x], [yB, evap_y+evap_h/2], color=pipe_c, linewidth=2.5)
ax.plot([xL, evap_x], [evap_y+evap_h/2, evap_y+evap_h/2], color=pipe_c, linewidth=2.5)
ax.plot([evap_x+evap_w, evap_x+evap_w], [yB, evap_y+evap_h/2], color=pipe_c, linewidth=2.5)
ax.plot([evap_x+evap_w, xR], [evap_y+evap_h/2, evap_y+evap_h/2], color=pipe_c, linewidth=2.5)

# Embedded directional arrows along loop: 1->2->3->4->1
arrows = [
    ((xR, 2.5), (xR, 5.2)),      # up right side (1->2)
    ((9.8, yT), (6.5, yT)),      # left across top (2->3)
    ((xL, 5.3), (xL, 2.5)),      # down left side (3->4)
    ((4.7, yB), (8.7, yB)),      # right across bottom (4->1)
]
for p0, p1 in arrows:
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=15, linewidth=2.2, color='black'))

# Labels for components
ax.text(comp_cx+1.0, comp_cy, 'Compressor', fontsize=12, va='center')
ax.text(cond_x + cond_w/2, cond_y + cond_h/2, 'Condenser', fontsize=12, ha='center', va='center', weight='bold')
ax.text(valve_x, valve_y-0.8, 'Expansion valve', fontsize=11, ha='center')
ax.text(evap_x + evap_w/2, evap_y + evap_h/2, 'Evaporator', fontsize=12, ha='center', va='center', weight='bold')

# Fan symbols near HX (simple 3-blade)
def fan(cx, cy, r=0.25, color='black'):
    ax.add_patch(Circle((cx, cy), r, facecolor='white', edgecolor=color, linewidth=1.2))
    for a in [20, 140, 260]:
        t = np.deg2rad(a)
        p1 = (cx + 0.02*np.cos(t), cy + 0.02*np.sin(t))
        p2 = (cx + 0.22*np.cos(t+0.5), cy + 0.22*np.sin(t+0.5))
        p3 = (cx + 0.12*np.cos(t-0.9), cy + 0.12*np.sin(t-0.9))
        ax.add_patch(Polygon([p1,p2,p3], closed=True, facecolor=color, edgecolor=color, linewidth=0.6))

fan(10.9, 6.35)
fan(10.9, 1.65)

# Heat arrows and ambient/cold space labels
ax.add_patch(FancyArrowPatch((10.9, 6.75), (10.9, 7.45), arrowstyle='-|>', mutation_scale=18, linewidth=2, color=heat_red))
ax.text(10.9, 7.58, 'Ambient, $T_0$', fontsize=11, ha='center')

ax.add_patch(FancyArrowPatch((10.9, 0.45), (10.9, 1.05), arrowstyle='-|>', mutation_scale=18, linewidth=2, color=heat_blue))
ax.text(10.9, 0.26, 'Cold space, $T_{CL}$', fontsize=11, ha='center')

# Compressor work arrow
ax.add_patch(FancyArrowPatch((12.8, 3.95), (12.15, 3.95), arrowstyle='-|>', mutation_scale=18, linewidth=2, color=gray))
ax.text(12.85, 4.25, '$W_{comp}$', fontsize=11, color=gray, ha='left')

# State markers
state_pos = {
    1: (10.45, 1.9),
    2: (10.45, 6.1),
    3: (2.95, 6.1),
    4: (2.95, 1.9),
}
phase_labels = {
    1: 'Low-P vapor',
    2: 'High-P vapor',
    3: 'High-P liquid',
    4: 'Low-P liquid + vapor',
}
for s, (x,y) in state_pos.items():
    ax.add_patch(Circle((x,y), 0.18, facecolor=state_fill, edgecolor='black', linewidth=1.2))
    ax.text(x, y, str(s), fontsize=10, ha='center', va='center', weight='bold')

ax.text(10.8, 1.92, phase_labels[1], fontsize=9, va='center')
ax.text(10.8, 6.12, phase_labels[2], fontsize=9, va='center')
ax.text(2.15, 6.12, phase_labels[3], fontsize=9, va='center', ha='right')
ax.text(2.15, 1.92, phase_labels[4], fontsize=9, va='center', ha='right')

# Equations
ax.text(12.0, 3.1, r'$W_{comp} = \dot{m}(h_2-h_1)$', fontsize=10, ha='left')
ax.text(7.0, 7.05, r'$Q_{cond} = \dot{m}(h_2-h_3)$', fontsize=10, ha='center')
ax.text(1.0, 3.3, r'$h_3=h_4$', fontsize=10, ha='left')
ax.text(7.0, 0.55, r'$Q_{evap} = \dot{m}(h_1-h_4)$', fontsize=10, ha='center')
ax.text(11.1, 0.9, r'$COP = \dfrac{Q_{evap}}{W_{comp}}$', fontsize=12, ha='left')

# Title
ax.text(7.0, 7.85, 'Single-Stage Vapor Compression Cycle Model', fontsize=18, ha='center', weight='bold')

out_png = '/Users/snarasi2/idaes-hvacr-cycles/vcrc_single_stage_schematic.png'
out_svg = '/Users/snarasi2/idaes-hvacr-cycles/vcrc_single_stage_schematic.svg'
fig.savefig(out_png, bbox_inches='tight', dpi=300)
fig.savefig(out_svg, bbox_inches='tight')
print(out_png)
print(out_svg)
