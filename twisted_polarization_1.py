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
cross_arrows = []
num_cross_arrows = 60

k_twist = 0.
k_wind = 0.
k_wave = 0.

omega_twist = 1.
omega_wind = 1.
omega_wave = 1.

phase_init_twist_pi = 0.
phase_init_wind_pi = 0.
phase_init_wave_pi = 0.

radius_wind = 0.
amp_wave_v = 1.
amp_wave_h = 1.

is_twist = False
is_wind = False
is_wave = False

is_reverse_twist = False
is_reverse_wind = False
is_wave_reverse = False

pos_neg_twist = 1.
pos_neg_wind = 1.
pos_neg_wave = 1.


""" Animation control """
is_play = False

""" Axis vectors """
vector_x_axis = np.array([1., 0., 0.])
vector_y_axis = np.array([0., 1., 0.])
vector_z_axis = np.array([0., 0., 1.])

""" Other parameters """


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
    def __init__(self, ax, origin, theta_twist, amplitude_h, amplitude_v):
        self.ax = ax
        self.origin = origin
        self.theta_twist = theta_twist
        self.amplitude_h = amplitude_h
        self.amplitude_v = amplitude_v

        # self.phase_h = phase_h
        # self.phase_v = phase_v

        self.vertical_axis = np.array([0., 0., 1.])
        self.horizontal_axis = np.array([0., 1., 0.])

        rot_matrix_x = Rotation.from_rotvec(self.theta_twist * vector_x_axis)
        self.vertical_axis_rotated = rot_matrix_x.apply(self.vertical_axis)
        self.horizontal_axis_rotated = rot_matrix_x.apply(self.horizontal_axis)

        x, y, z = self.origin[0], self.origin[1], self.origin[2]

        self.point_h = self.horizontal_axis * self.amplitude_h
        u, v, w = self.point_h[0], self.point_h[1], self.point_h[2]
        self.arrow_h = Arrow3d(self.ax, x, y, z, u, v, w, "green", 0.3, "-", 0.1)

        self.point_v = self.vertical_axis_rotated * self.amplitude_v
        u, v, w = self.point_v[0], self.point_v[1], self.point_v[2]
        self.arrow_v = Arrow3d(self.ax, x, y, z, u, v, w, "red", 1, "-", 0.1)

    def update_diagrams(self):
        rot_matrix_x = Rotation.from_rotvec(self.theta_twist * vector_x_axis)
        self.vertical_axis_rotated = rot_matrix_x.apply(self.vertical_axis)
        self.horizontal_axis_rotated = rot_matrix_x.apply(self.horizontal_axis)

        self.point_h = self.horizontal_axis_rotated * self.amplitude_h
        self.point_v = self.vertical_axis_rotated * self.amplitude_v

        x, y, z = self.origin[0], self.origin[1], self.origin[2]

        self.arrow_h.set_origin(x, y, z)
        u, v, w = self.point_h[0], self.point_h[1], self.point_h[2]
        self.arrow_h.set_vector(u, v, w)

        self.arrow_v.set_origin(x, y, z)
        u, v, w = self.point_v[0], self.point_v[1], self.point_v[2]
        self.arrow_v.set_vector(u, v, w)

    def set_amplitude(self, value_h, value_v):
        self.amplitude_h, self.amplitude_v = value_h, value_v
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


def set_is_twist(value):
    global is_twist
    is_twist = value
    update_diagrams()


def set_turn_twist(value):
    global k_twist
    k_twist = value
    update_diagrams()


def set_phase_init_twist(value):
    global phase_init_twist_pi
    phase_init_twist_pi = value
    update_diagrams()


def set_reverse_twist(value):
    global pos_neg_twist
    if value:
        pos_neg_twist = -1.
    else:
        pos_neg_twist = 1.


def set_is_wind(value):
    global is_wind
    is_wind = value
    if not is_wind:
        for i in range(num_cross_arrows):
            x = i * (x_max - x_min) / num_cross_arrows
            wind_h = 0.
            wind_v = 0.
            cross_arrows[i].set_origin(x, wind_h, wind_v)

        wind_h = x_seq * 0.
        wind_v = x_seq * 0.
        plt_wind.set_data_3d(x_seq, wind_h, wind_v)

    update_diagrams()


def set_radius_wind(value):
    global radius_wind
    radius_wind = value
    update_diagrams()


def set_turn_wind(value):
    global k_wind
    k_wind = value
    update_diagrams()


def set_phase_init_wind(value):
    global phase_init_wind_pi
    phase_init_wind_pi = value
    update_diagrams()


def set_reverse_wind(value):
    global pos_neg_wind
    if value:
        pos_neg_wind = -1.
    else:
        pos_neg_wind = 1.


def set_is_wave(value):
    global is_wave
    is_wave = value
    update_diagrams()


def set_k_wave(value):
    global k_wave
    k_wave = value
    update_diagrams()


def set_phase_init_wave(value):
    global phase_init_wave_pi
    phase_init_wave_pi = value
    update_diagrams()


def set_reverse_wave(value):
    global pos_neg_wave
    if value:
        pos_neg_wave = -1.
    else:
        pos_neg_wave = 1.


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


def clear_snapped_path():
    path_v_snap.clear_path()
    path_h_snap.clear_path()
    path_helicity_snap.clear_path()


def create_parameter_setter():
    # Twist
    frm_twist = ttk.Labelframe(root, relief='ridge', text="Twist", labelanchor='n', width=100)
    frm_twist.pack(side='left')

    var_is_twist = tk.IntVar(root)
    chk_is_twist = tk.Checkbutton(frm_twist, text="", variable=var_is_twist,
                                  command=lambda: set_is_twist(var_is_twist.get()))
    chk_is_twist.pack(side='left')
    var_is_twist.set(is_twist)

    lbl_turn_twist = tk.Label(frm_twist, text="Turn")
    lbl_turn_twist.pack(side='left')
    var_turn_twist = tk.StringVar(root)
    var_turn_twist.set(str(k_twist))
    spn_turn_twist = tk.Spinbox(
        frm_twist, textvariable=var_turn_twist, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_turn_twist(float(var_turn_twist.get())), width=8
    )
    spn_turn_twist.pack(side='left')

    lbl_phase_init_twist = tk.Label(frm_twist, text="Ph.")
    lbl_phase_init_twist.pack(side='left')
    var_phase_init_twist = tk.StringVar(root)
    var_phase_init_twist.set(str(phase_init_twist_pi))
    spn_phase_init_twist = tk.Spinbox(
        frm_twist, textvariable=var_phase_init_twist, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_phase_init_twist(float(var_phase_init_twist.get())), width=8
    )
    spn_phase_init_twist.pack(side='left')

    var_is_reverse_twist = tk.IntVar(root)
    chk_is_reverse_twist = tk.Checkbutton(frm_twist, text="Reverse", variable=var_is_reverse_twist,
                                          command=lambda: set_reverse_twist(var_is_reverse_twist.get()))
    chk_is_reverse_twist.pack(side='left')
    var_is_reverse_twist.set(is_reverse_twist)

    # Wind
    frm_wind = ttk.Labelframe(root, relief='ridge', text="Wind", labelanchor='n', width=100)
    frm_wind.pack(side='left')

    var_is_wind = tk.IntVar(root)
    chk_is_wind = tk.Checkbutton(frm_wind, text="", variable=var_is_wind,
                                 command=lambda: set_is_wind(var_is_wind.get()))
    chk_is_wind.pack(side='left')
    var_is_wind.set(is_wind)

    lbl_radius_wind = tk.Label(frm_wind, text="R")
    lbl_radius_wind.pack(side='left')
    var_radius_wind = tk.StringVar(root)
    var_radius_wind.set(str(radius_wind))
    spn_radius_wind = tk.Spinbox(
        frm_wind, textvariable=var_radius_wind, format='%.1f', from_=0., to=4., increment=0.1,
        command=lambda: set_radius_wind(float(var_radius_wind.get())), width=8
    )
    spn_radius_wind.pack(side='left')

    lbl_turn_wind = tk.Label(frm_wind, text="Turn")
    lbl_turn_wind.pack(side='left')
    var_turn_wind = tk.StringVar(root)
    var_turn_wind.set(str(k_wind))
    spn_turn_wind = tk.Spinbox(
        frm_wind, textvariable=var_turn_wind, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_turn_wind(float(var_turn_wind.get())), width=8
    )
    spn_turn_wind.pack(side='left')

    lbl_phase_init_wind = tk.Label(frm_wind, text="Ph.")
    lbl_phase_init_wind.pack(side='left')
    var_phase_init_wind = tk.StringVar(root)
    var_phase_init_wind.set(str(phase_init_twist_pi))
    spn_phase_init_wind = tk.Spinbox(
        frm_wind, textvariable=var_phase_init_wind, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_phase_init_wind(float(var_phase_init_wind.get())), width=8
    )
    spn_phase_init_wind.pack(side='left')

    var_is_reverse_wind = tk.IntVar(root)
    chk_is_reverse_wind = tk.Checkbutton(frm_wind, text="Reverse", variable=var_is_reverse_wind,
                                         command=lambda: set_reverse_wind(var_is_reverse_wind.get()))
    chk_is_reverse_wind.pack(side='left')
    var_is_reverse_wind.set(is_reverse_wind)

    # Wave
    frm_wave = ttk.Labelframe(root, relief='ridge', text="Wave", labelanchor='n', width=100)
    frm_wave.pack(side='left')

    var_is_wave = tk.IntVar(root)
    chk_is_wave = tk.Checkbutton(frm_wave, text="", variable=var_is_wave,
                                 command=lambda: set_is_wave(var_is_wave.get()))
    chk_is_wave.pack(side='left')
    var_is_wave.set(is_wave)

    lbl_k_wave = tk.Label(frm_wave, text="k")
    lbl_k_wave.pack(side='left')
    var_k_wave = tk.StringVar(root)
    var_k_wave.set(str(k_wave))
    spn_k_wave = tk.Spinbox(
        frm_wave, textvariable=var_k_wave, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_k_wave(float(var_k_wave.get())), width=8
    )
    spn_k_wave.pack(side='left')

    lbl_phase_init_wave = tk.Label(frm_wave, text="Phase")
    lbl_phase_init_wave.pack(side='left')
    var_phase_init_wave = tk.StringVar(root)
    var_phase_init_wave.set(str(phase_init_wave_pi))
    spn_phase_init_wave = tk.Spinbox(
        frm_wave, textvariable=var_phase_init_wave, format='%.1f', from_=-4., to=4., increment=0.1,
        command=lambda: set_phase_init_wave(float(var_phase_init_wave.get())), width=8
    )
    spn_phase_init_wave.pack(side='left')

    var_is_reverse_wave = tk.IntVar(root)
    chk_is_reverse_wave = tk.Checkbutton(frm_wave, text="Reverse", variable=var_is_reverse_wave,
                                         command=lambda: set_reverse_wave(var_is_reverse_wave.get()))
    chk_is_reverse_wave.pack(side='left')
    var_is_reverse_wave.set(is_reverse_wind)

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


def set_path():
    path_h.clear_path()
    path_v.clear_path()
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


def reset():
    global is_play
    cnt.reset()
    update_diagrams()


def switch():
    global is_play
    is_play = not is_play


def update_diagrams():
    global wind_h, wind_v
    t = cnt.get() / 10.

    for i in range(num_cross_arrows):
        x = i * (x_max - x_min) / num_cross_arrows

        if is_twist:
            phase_twist = (k_twist * x - pos_neg_twist * omega_twist * t + phase_init_twist_pi) * np.pi
            cross_arrows[i].set_twist(phase_twist)

        if is_wind:
            phase_wind = (k_wind * x - omega_wind * t + phase_init_wind_pi) * np.pi
            wind_h = pos_neg_wind * radius_wind * np.cos(phase_wind)
            wind_v = radius_wind * np.sin(phase_wind)
            cross_arrows[i].set_origin(x, wind_h, wind_v)

            wind_h = pos_neg_wind * radius_wind * np.cos((k_wind * x_seq - omega_wind * t + phase_init_wind_pi) * np.pi)
            wind_v = radius_wind * np.sin((k_wind * x_seq - omega_wind * t + phase_init_wind_pi) * np.pi)
            plt_wind.set_data_3d(x_seq, wind_h, wind_v)

        if is_wave:
            phase_wave = (k_wave * x - omega_wave * t + phase_init_wave_pi) * np.pi
            amp_h = pos_neg_wave * np.cos(phase_wave)
            amp_v = np.sin(phase_wave)
            cross_arrows[i].set_amplitude(amp_h, amp_v)

    set_path()


def update(f):
    if is_play:
        cnt.count_up()
        update_diagrams()


""" main loop """
if __name__ == "__main__":
    cnt = Counter(ax=ax0, is3d=True, xy=np.array([x_min, y_max]), z=z_max, label="Step=")
    draw_static_diagrams()
    create_animation_control()
    create_parameter_setter()

    t = 0.
    x_seq = np.arange(x_min, x_max, 0.005)

    for i_ in range(num_cross_arrows):
        x_ = i_ * (x_max - x_min) / num_cross_arrows
        origin = np.array([x_, 0., 0.])

        cross_arrow = CrossArrow(ax0, origin, 0, amp_wave_h, amp_wave_v)

        phase_twist = (k_twist * x_ - omega_twist * t + phase_init_twist_pi) * np.pi
        cross_arrow.set_twist(phase_twist)

        phase_wind = (k_wind * x_ - omega_wind * t + phase_init_wind_pi) * np.pi
        wind_h = radius_wind * np.cos(phase_wind)
        wind_v = radius_wind * np.sin(phase_wind)
        cross_arrow.set_origin(x_, wind_h, wind_v)

        if is_wave:
            phase_wave = (k_wave * x_ - omega_wave * t + phase_init_wave_pi) * np.pi
            amp_h = np.cos(phase_wave)
            amp_v = np.sin(phase_wave)
            cross_arrow.set_amplitude(amp_h, amp_v)

        cross_arrows.append(cross_arrow)

    path_h = Path(ax0, "-", 0.5, "green")
    path_v = Path(ax0, "-", 1, "red")
    path_v_plus_h = Path(ax0, "-", 2, "orange")
    set_path()

    wind_h = radius_wind * np.cos((k_wind * x_seq - omega_wind * t + phase_init_wind_pi) * np.pi)
    wind_v = radius_wind * np.sin((k_wind * x_seq - omega_wind * t + phase_init_wind_pi) * np.pi)

    plt_wind, = ax0.plot(x_seq, wind_h, wind_v, color="gray", ls="--", linewidth=1)

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
