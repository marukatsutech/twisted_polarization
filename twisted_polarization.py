""" Twisted polarization """
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import proj3d
from scipy.spatial.transform import Rotation

""" Global variables """

""" Animation control """
is_play = False

""" Axis vectors """
vector_x_axis = np.array([1., 0., 0.])
vector_y_axis = np.array([0., 1., 0.])
vector_z_axis = np.array([0., 0., 1.])

""" Other parameters """
num_turn = 0.

t = 0.
omega = 1.
r_twist = 0.
phase_twist_pi = 0.

cross_arrows = []
num_cross_arrows = 80

amplitude_v = 1.
amplitude_h = 1.

is_reverse = False
pos_neg = 1.

is_wave = True

""" Create figure and axes """
title_ax0 = "Twisted polarization"
# title_ax1 = "aaa"
title_tk = title_ax0

x_min = 0.
x_max = 4.
y_min = -2.
y_max = 2.
z_min = -2.
z_max = 2.

fig = Figure()
# fig = Figure(facecolor='black')
ax0 = fig.add_subplot(111, projection='3d')
ax0.set_box_aspect((2, 1, 1))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel(r"$z(\pi)$")
ax0.set_ylabel("x")
ax0.set_zlabel("y")
ax0.set_xlim(x_min, x_max)
ax0.set_ylim(y_min, y_max)
ax0.set_zlim(z_min, z_max)

# ax0.set_facecolor('black')
# ax0.axis('off')

""" Embed in Tkinter """
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

""" Global objects of Tkinter """
var_snap_op = tk.IntVar()

""" Classes and functions """


class Counter:
    def __init__(self, is3d=None, ax=None, xy=None, z=None, label="", color=None):
        self.is3d = is3d if is3d is not None else False
        self.ax = ax
        self.x, self.y = xy[0], xy[1]
        self.z = z if z is not None else 0
        self.label = label
        self.color = color

        self.count = 0

        if not is3d:
            self.txt_step = self.ax.text(self.x, self.y, self.label + str(self.count), color=color)
        else:
            self.txt_step = self.ax.text2D(self.x, self.y, self.label + str(self.count), color=color)
            self.xz, self.yz, _ = proj3d.proj_transform(self.x, self.y, self.z, self.ax.get_proj())
            self.txt_step.set_position((self.xz, self.yz))

    def count_up(self):
        self.count += 1
        self.txt_step.set_text(self.label + str(self.count))

    def reset(self):
        self.count = 0
        self.txt_step.set_text(self.label + str(self.count))

    def get(self):
        return self.count


class Arrow3d:
    def __init__(self, ax, x, y, z, u, v, w, color, line_width, line_style, arrow_length_ratio):
        self.ax = ax
        self.x, self.y, self.z = x, y, z
        self.u, self.v, self.w = u, v, w
        self.color = color
        self.line_width = line_width
        self.line_style = line_style
        self.arrow_length_ratio = arrow_length_ratio

        self.qvr = self.ax.quiver(self.x, self.y, self.z, self.u, self.v, self.w,
                                  length=1, color=self.color, normalize=False,
                                  linewidth=self.line_width, linestyle=self.line_style,
                                  arrow_length_ratio=self.arrow_length_ratio)

    def _update_quiver(self):
        self.qvr.remove()
        self.qvr = self.ax.quiver(self.x, self.y, self.z, self.u, self.v, self.w,
                                  length=1, color=self.color, normalize=False,
                                  linewidth=self.line_width, linestyle=self.line_style,
                                  arrow_length_ratio=self.arrow_length_ratio)

    def set_origin(self, x, y, z):
        self.x, self.y, self.z = x, y, z
        self._update_quiver()

    def set_vector(self, u, v, w):
        self.u, self.v, self.w = u, v, w
        self._update_quiver()

    def get_origin(self):
        return np.array([self.x, self.y, self.z])

    def get_vector(self):
        return np.array([self.u, self.v, self.w])


class CrossArrow:
    def __init__(self, ax, origin, theta_twist, amplitude_v, amplitude_h):
        self.ax = ax
        self.origin = origin
        self.theta_twist = theta_twist
        self.amplitude_v = amplitude_v
        self.amplitude_h = amplitude_h
        # self.phase_v = phase_v
        # self.phase_h = phase_h

        self.vertical_axis = np.array([0., 0., 1.])
        self.horizontal_axis = np.array([0., 1., 0.])

        rot_matrix_x = Rotation.from_rotvec(self.theta_twist * vector_x_axis)
        self.vertical_axis_rotated = rot_matrix_x.apply(self.vertical_axis)
        self.horizontal_axis_rotated = rot_matrix_x.apply(self.horizontal_axis)

        x, y, z = self.origin[0], self.origin[1], self.origin[2]

        self.point_v = self.vertical_axis_rotated * self.amplitude_v
        u, v, w = self.point_v[0], self.point_v[1], self.point_v[2]
        self.arrow_v = Arrow3d(self.ax, x, y, z, u, v, w, "red", 1, "-", 0.1)

        self.point_h = self.horizontal_axis * self.amplitude_h
        u, v, w = self.point_h[0], self.point_h[1], self.point_h[2]
        self.arrow_h = Arrow3d(self.ax, x, y, z, u, v, w, "green", 0.3, "-", 0.1)

    def update_diagrams(self):
        rot_matrix_x = Rotation.from_rotvec(self.theta_twist * vector_x_axis)
        self.vertical_axis_rotated = rot_matrix_x.apply(self.vertical_axis)
        self.horizontal_axis_rotated = rot_matrix_x.apply(self.horizontal_axis)

        self.point_v = self.vertical_axis_rotated * self.amplitude_v
        self.point_h = self.horizontal_axis_rotated * self.amplitude_h

        x, y, z = self.origin[0], self.origin[1], self.origin[2]

        self.arrow_v.set_origin(x, y, z)
        u, v, w = self.point_v[0], self.point_v[1], self.point_v[2]
        self.arrow_v.set_vector(u, v, w)

        self.arrow_h.set_origin(x, y, z)
        u, v, w = self.point_h[0], self.point_h[1], self.point_h[2]
        self.arrow_h.set_vector(u, v, w)

    def set_amplitude(self, value_v, value_h):
        self.amplitude_v, self.amplitude_h = value_v, value_h
        self.update_diagrams()

    def set_origin(self, x, y, z):
        self.origin[0], self.origin[1], self.origin[2] = x, y, z
        self.update_diagrams()

    def set_twist(self, value):
        self.theta_twist = value
        self.update_diagrams()

    def get_vector_v(self):
        x, y, z = self.origin[0], self.origin[1], self.origin[2]
        u, v, w = self.point_v[0], self.point_v[1], self.point_v[2]
        return u + x, v + y, w + z

    def get_vector_h(self):
        x, y, z = self.origin[0], self.origin[1], self.origin[2]
        u, v, w = self.point_h[0], self.point_h[1], self.point_h[2]
        return u + x, v + y, w + z


class Path:
    def __init__(self, ax, line_style, line_width, color):
        self.ax = ax
        self.line_style = line_style
        self.line_width = line_width
        self.color = color

        self.is_draw_path = True

        self.x_path = []
        self.y_path = []
        self.z_path = []
        self.path, = self.ax.plot(np.array(self.x_path), np.array(self.y_path), np.array(self.z_path),
                                  color=self.color, linewidth=self.line_width, linestyle=self.line_style)

    def append_path(self, position):
        if self.is_draw_path:
            self.x_path.append(position[0])
            self.y_path.append(position[1])
            self.z_path.append(position[2])
            self.update_path()

    def update_path(self):
        self.path.set_data_3d(np.array(self.x_path), np.array(self.y_path), np.array(self.z_path))

    def clear_path(self):
        self.x_path = []
        self.y_path = []
        self.z_path = []
        self.update_path()

    def set_is_draw_path(self, value):
        self.is_draw_path = value


def set_path():
    path_v.clear_path()
    path_h.clear_path()
    path_v_plus_h.clear_path()
    for arrow in cross_arrows:
        p_v = arrow.get_vector_v()
        path_v.append_path(np.array([p_v[0], p_v[1], p_v[2]]))
        p_h = arrow.get_vector_h()
        path_h.append_path(np.array([p_h[0], p_h[1], p_h[2]]))

        path_v_plus_h.append_path(np.array([p_v[0],
                                            1. / np.sqrt(2.) * (p_v[1] + p_h[1]),
                                            1. / np.sqrt(2.) * (p_v[2] + p_h[2])
                                            ]))


def set_path_snap():
    path_v_snap.clear_path()
    path_h_snap.clear_path()
    path_helicity_snap.clear_path()
    for arrow in cross_arrows:
        p_v = arrow.get_vector_v()
        path_v_snap.append_path(np.array([p_v[0], p_v[1], p_v[2]]))
        p_h = arrow.get_vector_h()
        path_h_snap.append_path(np.array([p_h[0], p_h[1], p_h[2]]))
        path_helicity_snap.append_path(np.array([p_v[0],
                                       1. / np.sqrt(2.) * (p_v[1] + p_h[1]),
                                       1. / np.sqrt(2.) * (p_v[2] + p_h[2])
                                       ]))


def set_turn(value):
    global num_turn
    num_turn = value
    update_diagrams()


def set_radius(value):
    global r_twist
    r_twist = value
    update_diagrams()


def set_phase(value):
    global phase_twist_pi
    phase_twist_pi = value
    update_diagrams()


def clear_snapped_path():
    path_v_snap.clear_path()
    path_h_snap.clear_path()


def set_reverse_twist(value):
    global pos_neg
    if value:
        pos_neg = -1.
    else:
        pos_neg = 1.


def set_is_wave(value):
    global is_wave
    is_wave = value
    update_diagrams()


def create_parameter_setter():
    frm_twist = ttk.Labelframe(root, relief='ridge', text="Twist", labelanchor='n', width=100)
    frm_twist.pack(side='left')

    lbl_radius = tk.Label(frm_twist, text="Radius")
    lbl_radius.pack(side='left')
    var_radius = tk.StringVar(root)
    var_radius.set(str(r_twist))
    spn_radius = tk.Spinbox(
        frm_twist, textvariable=var_radius, format='%.1f', from_=0., to=4., increment=0.1,
        command=lambda: set_radius(float(var_radius.get())), width=8
    )

    spn_radius.pack(side='left')
    lbl_turn = tk.Label(frm_twist, text="Turn(2pi/turn)")
    lbl_turn.pack(side='left')
    var_turn = tk.StringVar(root)
    var_turn.set(str(num_turn))
    spn_turn = tk.Spinbox(
        frm_twist, textvariable=var_turn, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_turn(float(var_turn.get())), width=8
    )
    spn_turn.pack(side='left')

    lbl_phase = tk.Label(frm_twist, text="Phase(pi)")
    lbl_phase.pack(side='left')
    var_phase = tk.StringVar(root)
    var_phase.set(str(phase_twist_pi))
    spn_phase = tk.Spinbox(
        frm_twist, textvariable=var_phase, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_phase(float(var_phase.get())), width=8
    )
    spn_phase.pack(side='left')

    var_is_reverse = tk.IntVar(root)
    chk_is_reverse = tk.Checkbutton(frm_twist, text="Reverse rotation", variable=var_is_reverse,
                                    command=lambda: set_reverse_twist(var_is_reverse.get()))
    chk_is_reverse.pack(side='left')
    var_is_reverse.set(is_reverse)

    frm_snap = ttk.Labelframe(root, relief='ridge', text="Snap path", labelanchor='n', width=100)
    frm_snap.pack(side='left')

    # var_snap_op = tk.IntVar()
    rd_op_snap = tk.Radiobutton(frm_snap, text="Snap", value=1, variable=var_snap_op,
                                command=lambda: set_path_snap())
    rd_op_snap.pack(side='left')

    rd_op_clear = tk.Radiobutton(frm_snap, text="Clear", value=2, variable=var_snap_op,
                                 command=lambda: clear_snapped_path())
    rd_op_clear.pack(side='left')
    var_snap_op.set(2)

    frm_wave = ttk.Labelframe(root, relief="ridge", text="Wave", labelanchor='n')
    frm_wave.pack(side="left", fill=tk.Y)

    var_chk_wave = tk.BooleanVar(root)
    chk_wave = tk.Checkbutton(frm_wave, text="Apply", variable=var_chk_wave,
                              command=lambda: set_is_wave(var_chk_wave.get()))
    chk_wave.pack(anchor=tk.W)
    var_chk_wave.set(is_wave)


def create_animation_control():
    frm_anim = ttk.Labelframe(root, relief='ridge', text="Animation", labelanchor='n')
    frm_anim.pack(side='left', fill=tk.Y)
    btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
    btn_play.pack(side='left')
    btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
    btn_reset.pack(side='left')
    # btn_clear = tk.Button(frm_anim, text="Clear path", command=lambda: aaa())
    # btn_clear.pack(fill=tk.X)


def create_center_lines(ax, x_min, x_max, y_min, y_max, z_min, z_max):
    line_axis_x = art3d.Line3D([x_min, x_max], [0., 0.], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax.add_line(line_axis_x)
    line_axis_y = art3d.Line3D([0., 0.], [y_min, y_max], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax.add_line(line_axis_y)
    line_axis_z = art3d.Line3D([0., 0.], [0., 0.], [z_min, z_max], color="gray", ls="-.", linewidth=1)
    ax.add_line(line_axis_z)


def draw_static_diagrams():
    create_center_lines(ax0, x_min, x_max, y_min, y_max, z_min, z_max)


def set_twist():
    for i in range(num_cross_arrows):
        theta_delta = np.pi * num_turn
        x = i * (x_max - x_min) / num_cross_arrows
        theta = x * theta_delta
        cross_arrows[i].set_twist(theta)


def update_diagrams():
    global t, rot_1_h, rot_1_v
    t = cnt.get() / 10.

    rot_1_v = pos_neg * r_twist * np.sin((x_seq - omega * t + phase_twist_pi) * np.pi * num_turn)
    rot_1_h = r_twist * np.cos((x_seq - omega * t + phase_twist_pi) * np.pi * num_turn)
    plt_rot_1.set_data_3d(x_seq, rot_1_h, rot_1_v)

    set_twist()

    for i in range(num_cross_arrows):
        x = i * (x_max - x_min) / num_cross_arrows
        if is_wave:
            amp_v = amplitude_v * np.sin((x - omega * t) * np.pi)
            amp_h = amplitude_h * np.cos((x - omega * t) * np.pi)
        else:
            amp_v, amp_h = 1., 1.
        cross_arrows[i].set_amplitude(amp_v, amp_h)
        index = np.argmin(np.abs(x - x_seq))
        z = rot_1_v[index]
        y = rot_1_h[index]
        cross_arrows[i].set_origin(x, y, z)

    set_path()


def reset():
    global is_play
    cnt.reset()
    update_diagrams()


def switch():
    global is_play
    is_play = not is_play


def update(f):
    parameter = fr"$Twist:{num_turn}/2\pi, Phase:{phase_twist_pi}(\pi)$"
    txt_parameter.set_text(parameter)
    if is_play:
        cnt.count_up()
        update_diagrams()


""" main loop """
if __name__ == "__main__":
    cnt = Counter(ax=ax0, is3d=True, xy=np.array([x_min, y_max]), z=z_max, label="Step=")
    draw_static_diagrams()
    create_animation_control()
    create_parameter_setter()

    x_seq = np.arange(x_min, x_max, 0.005)

    rot_1_v = r_twist * np.sin((x_seq - omega * t + phase_twist_pi) * np.pi * num_turn)
    rot_1_h = r_twist * np.cos((x_seq - omega * t + phase_twist_pi) * np.pi * num_turn)

    plt_rot_1, = ax0.plot(x_seq, rot_1_h, rot_1_v, color="gray", ls="--", linewidth=1)

    for i_ in range(num_cross_arrows):
        x_ = i_ * (x_max - x_min) / num_cross_arrows
        origin = np.array([x_, 0., 0.])
        amp_v = amplitude_v * np.sin((x_ - omega * t) * np.pi)
        amp_h = amplitude_h * np.cos((x_ - omega * t) * np.pi)

        cross_arrow = CrossArrow(ax0, origin, 0, amp_v, amp_h)
        cross_arrows.append(cross_arrow)

    path_v = Path(ax0, "-", 1, "red")
    path_h = Path(ax0, "-", 0.5, "green")
    path_v_plus_h = Path(ax0, "-", 2, "orange")
    set_path()

    parameter = fr"$Twist:{num_turn}/2\pi, Phase:{phase_twist_pi}(\pi)$"

    txt_parameter = ax0.text2D(x_min, y_max, parameter, fontsize="12")
    xz, yz, _ = proj3d.proj_transform(x_min + 0.5, y_max, z_max - 0.3, ax0.get_proj())
    txt_parameter.set_position((xz, yz))

    dummy1, = ax0.plot(np.array([0, 0]), np.array([0, 0]), np.array([0, 0]),
                       color="red", linewidth=1, linestyle="-", label="Vertical component(V)")
    dummy2, = ax0.plot(np.array([0, 0]), np.array([0, 0]), np.array([0, 0]),
                       color="green", linewidth=0.3, linestyle="-", label="Horizontal component(H)")
    dummy3, = ax0.plot(np.array([0, 0]), np.array([0, 0]), np.array([0, 0]),
                       color="orange", linewidth=2, linestyle="-", label=r"$Helicity(\frac{1}{\sqrt{2}}(V + H))$")

    path_v_snap = Path(ax0, "--", 1, "red")
    path_h_snap = Path(ax0, "--", 0.5, "green")
    path_helicity_snap = Path(ax0, "--", 1, "orange")

    ax0.legend(loc='lower right', fontsize=8)
    # ax1.legend(loc='lower right', fontsize=8)

    anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
    root.mainloop()
