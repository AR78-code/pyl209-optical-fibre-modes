"""Interactive Matplotlib demo; also usable with %matplotlib widget in Jupyter."""
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).parent / ".mplconfig"))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.patches import Circle
from fibre import Simulation


def create_demo():
    sim = Simulation()
    r = np.linspace(0, 12.5, 600)
    x = np.linspace(-12.5,12.5,241)
    xx,yy = np.meshgrid(x,x)
    rr = np.hypot(xx,yy)
    pp = np.arctan2(yy,xx)
    phi = np.linspace(0,2*np.pi,241)
    fig, ax = plt.subplots(2,3,figsize=(15,9))
    fig.subplots_adjust(bottom=.19,top=.88,wspace=.36,hspace=.42)
    fig.suptitle("Assignment 5 · Near-field fibre modes\nGhatak §§13.5–13.6 · coherent centred Gaussian launch",fontsize=16)
    extent=[-12.5,12.5,-12.5,12.5]
    im_out=ax[0,0].imshow(np.zeros_like(rr),origin="lower",extent=extent,cmap="magma",vmin=0,vmax=.25)
    ax[0,0].add_patch(Circle((0,0),sim.fibre.a,fill=False,ec="cyan",lw=.8))
    ax[0,0].set(title="Output pattern |E(x,y;L)|²",xlabel="x (µm)",ylabel="y (µm)")
    fig.colorbar(im_out,ax=ax[0,0],label="Intensity / lens-passed power (µm⁻²)")
    im_polar=ax[0,1].imshow(np.zeros((241,len(r))),origin="lower",aspect="auto",
                         extent=[0,12.5,0,2*np.pi],cmap="magma",vmin=0,vmax=.25)
    ax[0,1].set(title="Output |E(r,φ;L)|²",xlabel="r (µm)",ylabel="φ (rad)")
    ax[0,1].set_yticks([0,np.pi,2*np.pi],["0","π","2π"])
    line,=ax[0,2].plot(r,r*0,label="Guided output")
    launch,=ax[0,2].plot(r,r*0,"--",label="Focused input")
    ax[0,2].set(title="Radial intensity profiles",xlabel="r (µm)",ylabel="Intensity / power (µm⁻²)")
    ax[0,2].legend(fontsize=8)
    mode_img=ax[1,0].imshow(np.zeros_like(rr),origin="lower",extent=extent,cmap="RdBu_r",vmin=-1,vmax=1)
    ax[1,0].add_patch(Circle((0,0),sim.fibre.a,fill=False,ec="k",lw=.7))
    ax[1,0].set(xlabel="x (µm)",ylabel="y (µm)")
    fig.colorbar(mode_img,ax=ax[1,0],label="Normalized modal field (µm⁻¹)")
    mode_polar=ax[1,1].imshow(np.zeros((len(phi),len(r))),origin="lower",aspect="auto",
                           extent=[0,12.5,0,2*np.pi],cmap="RdBu_r",vmin=-1,vmax=1)
    ax[1,1].set(title="Modal field ψ(r,φ)",xlabel="r (µm)",ylabel="φ (rad)")
    ax[1,1].set_yticks([0,np.pi,2*np.pi],["0","π","2π"])
    bars=ax[1,2].bar(np.arange(len(sim.active)),np.zeros(len(sim.active)))
    ax[1,2].set(xticks=np.arange(len(sim.active)),xticklabels=[m.label for m in sim.active],
                title="Excited modal power fractions",ylabel="Fraction of lens-passed power",ylim=(0,1))
    ax[1,2].tick_params(axis="x",labelsize=7)
    length=Slider(fig.add_axes([.15,.11,.65,.023]),"L (mm)",0,5,valinit=.25)
    ratio=Slider(fig.add_axes([.15,.07,.65,.023]),"D/f",.1,1,valinit=1)
    mode=Slider(fig.add_axes([.15,.03,.65,.023]),"Mode index",0,len(sim.modes)-1,valinit=0,valstep=1)
    status=fig.text(.5,.915,"",ha="center",fontsize=10)
    saved=fig.add_axes([.84,.065,.12,.045]); button=Button(saved,"Save snapshot")
    cache={}
    def update(_=None):
        q=ratio.val
        if cache.get("ratio")!=q:
            cache.update(ratio=q,c=sim.coefficients(q),input=sim.focused_field(r,q))
        c=cache["c"]
        out=abs(sim.output(r,length.val,q,c))**2
        im_out.set_data(abs(sim.output(rr,length.val,q,c))**2)
        im_polar.set_data(np.tile(out,(len(phi),1)))
        line.set_ydata(out); launch.set_ydata(abs(cache["input"])**2)
        ax[0,2].relim(); ax[0,2].autoscale_view()
        for bar,weight in zip(bars,abs(c)**2): bar.set_height(weight)
        selected=sim.modes[int(mode.val)]
        v=selected.field(rr,pp); maxv=np.max(abs(v))
        mode_img.set_data(v); mode_img.set_clim(-maxv,maxv)
        mode_polar.set_data(selected.field(r[None,:],phi[:,None]))
        mode_polar.set_clim(-maxv,maxv)
        ax[1,0].set_title(f"{selected.label} · cosine orientation")
        status.set_text(f"V = {sim.fibre.V:.4f}  |  Guided fraction = {sum(abs(c)**2):.4%}  |  "
                        f"f = {20/q:.2f} mm  |  Selected n_eff = {selected.neff:.7f}")
        fig.canvas.draw_idle()
    for slider in [length,ratio,mode]: slider.on_changed(update)
    def save(_):
        Path("results").mkdir(exist_ok=True)
        fig.savefig("results/interactive_snapshot.png",dpi=170)
    button.on_clicked(save)
    update()
    # Keep widget references alive after returning from this function.
    fig._demo = dict(sim=sim,length=length,ratio=ratio,mode=mode,button=button,update=update)
    return fig


if __name__ == "__main__":
    create_demo()
    plt.show()
