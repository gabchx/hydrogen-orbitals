import math  # For math.factorial instead of np.math.factorial
import numpy as np
from scipy.special import eval_genlaguerre, sph_harm
import matplotlib.pyplot as plt
import matplotlib.animation as animation

a = 0.529e-10

def radial(r, n, l):
    p = 2*r / (n*a)
    lag_gen = eval_genlaguerre(n-l-1, 2*l+1, p)
    # Use math.factorial instead of np.math.factorial
    return (
        np.sqrt((2/(n*a))**3
                * math.factorial(n-l-1)
                / (2*n*math.factorial(n+l)))
        * np.exp(-p/2) * p**l * lag_gen
    )

def wave_func(r, theta, phi, nlm):
    """
    r, theta, phi in spherical coords.
    sph_harm signature is sph_harm(m, l, phi, theta).
    """
    n, l, m = nlm
    return radial(r, n, l) * sph_harm(m, l, phi, theta)

def wave_func_cart(x, y, z, nlm):
    """
    Convert (x,y,z) -> (r,theta,phi) and apply wave_func.
    """
    xy = x**2 + y**2
    r = np.sqrt(xy + z**2)
    theta = np.arctan2(np.sqrt(xy), z)
    phi = np.arctan2(y, x)
    return wave_func(r, theta, phi, nlm)


class Orbital_3D:
    def __init__(self, n, l, m, x=np.linspace(-2e-9, 2e-9, 201)):
        self.n, self.l, self.m = n, l, m
        self.x = x

        self.xyz, self.yxz, self.zxy = np.meshgrid(x, x, x, indexing='ij')
        psi_full = wave_func_cart(self.xyz, self.yxz, self.zxy, (n, l, m))
        prob_full = psi_full.real**2 + psi_full.imag**2

        norm = prob_full.sum()
        self.psi = psi_full / np.sqrt(norm)
        self.prob = prob_full / norm

class Transition_3D:
    """
    Stores a list of 'events' (either 'transition' or 'wait'), and
    then computes frames on-the-fly in save().
    """
    def __init__(self, orb, fps=30):
        self.orb = orb
        self.fps = fps
        self.events = []
        self.tot_frames = 0
        self.orbs = [(orb.n, orb.l, orb.m)]

    def transition(self, other, duration=1.0):
        """
        Record a transition from self.orb to 'other' over 'duration' seconds.
        """
        if isinstance(other, tuple):
            orb2 = Orbital_3D(*other)
        else:
            orb2 = other

        frames = int(duration * self.fps)
        self.events.append({
            'type': 'transition',
            'start': self.tot_frames,
            'end': self.tot_frames + frames,
            'orb1': self.orb,
            'orb2': orb2
        })
        self.orb = orb2
        self.orbs.append((orb2.n, orb2.l, orb2.m))
        self.tot_frames += frames

    def wait(self, duration):
        """
        Record a wait event (no change in wavefunction).
        """
        frames = int(duration * self.fps)
        self.events.append({
            'type': 'wait',
            'start': self.tot_frames,
            'end': self.tot_frames + frames,
            'orb': self.orb
        })
        self.tot_frames += frames

    def save(self, outfile="transition.mp4", samples=3000):
        """
        Render the animation by computing each frame on-the-fly.
        """
        orb = self.orb
        mid = orb.x.shape[0] // 2

        fig = plt.figure(figsize=(9,6))
        gs = fig.add_gridspec(2, 3, left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
        ax_3d = fig.add_subplot(gs[:, :2], projection='3d')
        axs = [fig.add_subplot(gs[k, 2]) for k in range(2)]

        ax_3d.view_init(22, -45)
        ax_3d.set_facecolor('black')

        # Flatten (x,y,z) for random sampling
        xi = orb.xyz.ravel()
        yi = orb.yxz.ravel()
        zi = orb.zxy.ravel()

        # Helper to lock axis but NOT set the title here
        def setup_3d(ax_3d):
            ax_3d.set_axis_off()
            ax_3d.set_xlim3d(orb.x[mid//2], orb.x[-mid//2])
            ax_3d.set_ylim3d(orb.x[mid//2], orb.x[-mid//2])
            ax_3d.set_zlim3d(orb.x[mid//2], orb.x[-mid//2])
            # Some lines for reference
            ax_3d.plot([orb.x[0],  orb.x[mid]],
                       [orb.x[-1], orb.x[-1]],
                       [orb.x[0],  orb.x[0]],
                       alpha=0.3, c="w", ls='--')
            ax_3d.plot([orb.x[0],  orb.x[0]],
                       [orb.x[mid], orb.x[-1]],
                       [orb.x[0],   orb.x[0]],
                       alpha=0.3, c="w", ls='--')
            ax_3d.plot([orb.x[0],  orb.x[0]],
                       [orb.x[-1], orb.x[-1]],
                       [orb.x[0],  orb.x[mid]],
                       alpha=0.3, c="w", ls='--')

        setup_3d(ax_3d)

        def find_event(frame_i):
            for e in self.events:
                if e['start'] <= frame_i < e['end']:
                    return e
            return self.events[-1]

        def animate(i):
            e = find_event(i)
            if e['type'] == 'wait':
                prob = e['orb'].prob
                # We show the label of the single orbital
                n, l, m = e['orb'].n, e['orb'].l, e['orb'].m
            else:
                # 'transition' event
                frames_total = e['end'] - e['start']
                frac = (i - e['start']) / frames_total
                angle = frac * (np.pi/2)
                c1, c2 = np.cos(angle), np.sin(angle)
                psi = c1*e['orb1'].psi + c2*e['orb2'].psi
                prob = psi.real**2 + psi.imag**2

                # We can always show the *target* orbital's label:
                # n, l, m = e['orb2'].n, e['orb2'].l, e['orb2'].m

                # Or if you want to "switch halfway" from old to new:
                if frac < 0.5:
                    n, l, m = e['orb1'].n, e['orb1'].l, e['orb1'].m
                else:
                    n, l, m = e['orb2'].n, e['orb2'].l, e['orb2'].m

            # We only re-scatter every other frame for performance
            if i % 2 == 0:
                ax_3d.clear()
                setup_3d(ax_3d)

                # Sample 3D points according to prob
                p = prob.ravel()
                p_sum = p.sum()
                if p_sum < 1e-16:
                    # fallback
                    p = np.ones_like(p)/p.size
                else:
                    p /= p_sum

                randices = np.random.choice(np.arange(p.size),
                                            samples, replace=True, p=p)
                x_rand, y_rand, z_rand = xi[randices], yi[randices], zi[randices]

                vmin, vmax = 0, prob.max()/5
                ax_3d.scatter(x_rand, y_rand, z_rand,
                              s=14,
                              c=prob.ravel()[randices],
                              cmap=plt.colormaps["cividis"],
                              alpha=0.25,
                              vmin=vmin, vmax=vmax)

            ax_3d.text2D(
                0.5, 0.95,                           # x,y position in the Axes coordinate [0..1]
                rf"$|n,l,m\rangle = |{n},{l},{m}\rangle$",
                transform=ax_3d.transAxes,
                ha="center", va="top",
                color="white", fontsize=16
            )


            # Update cross-section subplots each frame
            mid_slice_xy = prob[:, :, mid].T
            mid_slice_xz = prob[:, mid, :].T
            vmin, vmax = 0, prob.max()/5
            levels = np.linspace(vmin, vmax, 81)

            axs[0].clear()
            axs[0].set_axis_off()
            axs[0].set_title("$xy$-plane", color="w", y=0.02)
            axs[0].contourf(orb.x, orb.x, mid_slice_xy,
                            levels=levels, cmap=plt.colormaps["cividis"])

            axs[1].clear()
            axs[1].set_axis_off()
            axs[1].set_title("$xz$/$yz$-plane", color="w", y=0.02)
            axs[1].contourf(orb.x, orb.x, mid_slice_xz,
                            levels=levels, cmap=plt.colormaps["cividis"])

        ani = animation.FuncAnimation(
            fig, animate, frames=self.tot_frames,
            interval=1000/self.fps, repeat=False, blit=False
        )
        writer = animation.FFMpegWriter(fps=self.fps, codec="libx264", bitrate=4000)
        ani.save(outfile, writer=writer, dpi=120)
        print("Animation saved:", outfile)

if __name__ == "__main__":
    orb = Orbital_3D(1, 0, 0)     # Initial orbital n,l,m=1,0,0
    tr_3d = Transition_3D(orb, fps=30)

    tr_3d.wait(1.2)
    tr_3d.transition((2, 1, 0), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((2, 1, 1), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((3, 1, 0), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((3, 1, 1), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((3, 2, 0), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((4, 2, 0), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((4, 3, 0), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((4, 3, 1), duration=0.8)
    tr_3d.wait(1.2)
    tr_3d.transition((4, 1, 0), duration=0.8)
    tr_3d.wait(1.2)

    # Finally save the 3D animation
    tr_3d.save()

    # If you want to test the 2D transitions, you can do something like:
    # orb2d = Orbital(1, 0, 0)
    # tr_2d = Transition(orb2d, fps=60)
    # tr_2d.wait(1)
    # tr_2d.transition((2, 1, 1), duration=1)
    # tr_2d.save()
