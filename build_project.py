"""Rebuild figures, numerical tables, offline demo, and a runnable report notebook."""
import csv
import json
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR",str(Path(__file__).parent/".mplconfig"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import PowerNorm
import nbformat
from fibre import Simulation
from validate import run_validation

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                     "figure.facecolor":"white","savefig.facecolor":"white"})


def plot_results(s):
    x=np.linspace(-12.5,12.5,301); xx,yy=np.meshgrid(x,x)
    rr=np.hypot(xx,yy); pp=np.arctan2(yy,xx)
    extent=[-12.5,12.5,-12.5,12.5]
    def circle(ax):
        ax.add_patch(Circle((0,0),s.fibre.a,fill=False,ec="white",lw=.7,ls="--"))
        ax.set(xlabel="x (µm)",ylabel="y (µm)")
    c=s.coefficients(1)
    bound=.25
    fig,axes=plt.subplots(2,3,figsize=(12,8),layout="constrained")
    for ax,L in zip(axes.flat,[0,.05,.10,.25,.5,1.]):
        im=ax.imshow(abs(s.output(rr,L,1,c))**2,origin="lower",extent=extent,
                     cmap="magma",norm=PowerNorm(.5,vmin=0,vmax=bound))
        circle(ax); ax.set_title(f"L = {L:.2f} mm")
    fig.colorbar(im,ax=axes.ravel().tolist(),label="Intensity / lens-passed power (µm⁻²)",shrink=.85)
    fig.suptitle("Coherent near-field patterns as L varies · D/f = 1\nSix centred-launch modes; common square-root colour scale",fontsize=15)
    fig.savefig(OUT/"length_patterns.png",dpi=170);plt.close(fig)

    fig,axes=plt.subplots(1,4,figsize=(14,4),layout="constrained")
    for ax,q in zip(axes,[.1,.3,.5,1]):
        c=s.coefficients(q)
        im=ax.imshow(abs(s.output(rr,.25,q,c))**2,origin="lower",extent=extent,cmap="magma",norm=PowerNorm(.5,vmin=0,vmax=.25))
        circle(ax);ax.set_title(f"D/f = {q:.2f}\nGuided fraction {sum(abs(c)**2):.2%}")
    fig.colorbar(im,ax=axes.tolist(),label="Intensity / power (µm⁻²)",shrink=.8)
    fig.suptitle("Lens-ratio sweep at L = 0.25 mm · fixed D and incident beam\nCommon square-root colour scale",fontsize=14)
    fig.savefig(OUT/"lens_patterns.png",dpi=170);plt.close(fig)

    chosen=[next(m for m in s.modes if (m.ell,m.m)==lm)
            for lm in [(0,1),(0,2),(1,1),(2,1)]]
    r=np.linspace(0,12.5,400); phi=np.linspace(0,2*np.pi,241)
    fig,axes=plt.subplots(3,4,figsize=(13,10),layout="constrained")
    for col,m in enumerate(chosen):
        v=m.field(rr,pp); vmax=np.max(abs(v))
        im=axes[0,col].imshow(v,origin="lower",extent=extent,cmap="RdBu_r",vmin=-vmax,vmax=vmax)
        circle(axes[0,col]);axes[0,col].set_title(m.label+" · signed field")
        fig.colorbar(im,ax=axes[0,col],label="µm⁻¹",shrink=.7)
        im=axes[1,col].imshow(m.field(r[None,:],phi[:,None]),origin="lower",aspect="auto",
                             extent=[0,12.5,0,2*np.pi],cmap="RdBu_r",vmin=-vmax,vmax=vmax)
        axes[1,col].set(xlabel="r (µm)",ylabel="φ (rad)",yticks=[0,np.pi,2*np.pi],yticklabels=["0","π","2π"])
        axes[2,col].plot(r,m.radial(r),color="#126d78")
        axes[2,col].axvline(s.fibre.a,color="gray",ls="--",lw=.8)
        axes[2,col].set(xlabel="r (µm)",ylabel="R(r) (µm⁻¹)",title=f"n_eff = {m.neff:.6f}")
    fig.suptitle("Modal field distributions and profiles · cosine angular orientation",fontsize=15)
    fig.savefig(OUT/"modal_distributions.png",dpi=170);plt.close(fig)

    fig,axes=plt.subplots(9,5,figsize=(13,23),layout="constrained")
    for ax,m in zip(axes.flat,s.modes):
        v=m.field(rr,pp);peak=np.max(abs(v))
        ax.imshow(v,origin="lower",extent=extent,cmap="RdBu_r",vmin=-peak,vmax=peak)
        ax.add_patch(Circle((0,0),s.fibre.a,fill=False,ec="gray",lw=.5))
        ax.set_title(f"{m.label} · n_eff={m.neff:.5f}",fontsize=9)
        ax.set_xticks([-12.5,0,12.5]);ax.set_yticks([-12.5,0,12.5]);ax.tick_params(labelsize=7)
    fig.suptitle("All 45 supported LP families · signed field patterns\nAxes in µm; each panel scaled to its own field maximum; cosine orientation",fontsize=16)
    fig.savefig(OUT/"all_modes.png",dpi=140);plt.close(fig)

    radial=np.linspace(0,12.5,400); lengths=np.linspace(0,1,900); c=s.coefficients(1)
    intensity=np.array([abs(s.output(radial,L,1,c))**2 for L in lengths])
    fig,ax=plt.subplots(figsize=(10,4),layout="constrained")
    im=ax.imshow(intensity.T,origin="lower",aspect="auto",extent=[0,1,0,12.5],cmap="magma",norm=PowerNorm(.5,vmin=0,vmax=.25))
    ax.axhline(s.fibre.a,color="cyan",lw=.7,ls="--")
    ax.set(xlabel="L (mm)",ylabel="r (µm)",title="Radial near-field evolution · D/f = 1\nCoherent six-mode superposition · square-root colour scale")
    fig.colorbar(im,ax=ax,label="Intensity / power (µm⁻²)")
    fig.savefig(OUT/"length_evolution.png",dpi=180);plt.close(fig)


def export_tables(s):
    records=[m.record() for m in s.modes]
    with (OUT/"modes.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)
    with (OUT/"coupling_sweep.csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["D_over_f","f_mm","guided_fraction"]+
                                                [m.label+"_power_fraction" for m in s.active])
        for q in np.linspace(.1,1,91):
            powers=abs(s.coefficients(q))**2;writer.writerow([q,20/q,sum(powers),*powers])
    r=np.linspace(0,20,1001)
    for q,L in [(1,.25),(.3,.25),(1,0)]:
        out=s.output(r,L,q);inp=s.focused_field(r,q)
        np.savetxt(OUT/f"near_field_q{q:g}_L{L:g}mm.csv",np.c_[r,out.real,out.imag,abs(out)**2,abs(inp)**2],
                   delimiter=",",header="r_um,Re_E_per_um,Im_E_per_um,output_intensity_per_um2,input_intensity_per_um2",comments="")


def export_dashboard(s,validation):
    r=np.linspace(0,18,1801)
    modes=[dict(m.record(),radial=np.round(m.radial(r),11).tolist()) for m in s.modes]
    active=[s.modes.index(m) for m in s.active]
    launch=[]; bound=0
    basis=np.array([m.radial(r) for m in s.active])
    for q in np.linspace(.1,1,91):
        c=s.coefficients(q)
        bound=max(bound,float(np.max(np.sum(abs(c[:,None])*abs(basis),axis=0)**2)))
        launch.append(dict(c=np.c_[c.real,c.imag].tolist(),input=np.round(abs(s.focused_field(r,q))**2,11).tolist()))
    data=dict(r=np.round(r,8).tolist(),modes=modes,active=active,launch=launch,
              summary=s.summary(),validation=validation,intensity_bound=bound*1.01,
              waist_at_q1=s.fibre.wavelength*s.D/(np.pi*s.fill*s.D/2))
    template=(ROOT/"dashboard.template.html").read_text()
    (ROOT/"Fibre_Demo.html").write_text(template.replace("__DATA__",json.dumps(data,separators=(",",":"))))
    # A separate machine-readable export permits independent reproduction.
    np.savez_compressed(OUT/"radial_mode_data.npz",r_um=r,
                        radials=np.array([m.radial(r) for m in s.modes]),
                        beta=np.array([m.beta for m in s.modes]))


def create_notebook():
    md=(ROOT/"DERIVATION.md").read_text()
    cells=[nbformat.v4.new_markdown_cell(md),nbformat.v4.new_markdown_cell(
        "## 8. Numerical results\nRun from the project directory with the installed Python dependencies. "
        "The following cells compute the spectrum, coupling and validation; figures are regenerated by `build_project.py`."),
        nbformat.v4.new_code_cell("import numpy as np\nimport matplotlib.pyplot as plt\nfrom fibre import Simulation\nfrom validate import run_validation\nsim = Simulation()\nprint(f'V = {sim.fibre.V:.6f}; {len(sim.modes)} LP families; {len(sim.active)} excited radial modes.')"),
        nbformat.v4.new_code_cell("from IPython.display import display, Markdown, Image\nrows = ['| Mode | U | W | n_eff | Power fraction at D/f=1 |', '|---|---:|---:|---:|---:|']\nc = sim.coefficients(1.0)\nfor m, cj in zip(sim.active, c):\n    rows.append(f'| {m.label} | {m.u:.8f} | {m.w:.8f} | {m.neff:.9f} | {abs(cj)**2:.8f} |')\ndisplay(Markdown('\\n'.join(rows)))"),
        nbformat.v4.new_code_cell("report = run_validation(sim)\nassert report['all_passed'], report\nprint('Numerical checks passed.')")]
    for title,image in [("Modal fields in (r, φ), patterns and radial profiles","modal_distributions.png"),
                        ("Variable waveguide length L","length_patterns.png"),
                        ("Variable lens ratio D/f","lens_patterns.png")]:
        cells.extend([nbformat.v4.new_markdown_cell("### "+title),
                      nbformat.v4.new_code_cell(f"display(Image(filename='results/{image}'))")])
    cells.append(nbformat.v4.new_markdown_cell("## 9. Demonstration\nOpen `Fibre_Demo.html` and vary L and D/f. The two views show the summed output and individual modal fields."))
    nb=nbformat.v4.new_notebook(cells=cells,metadata={"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"}})
    nbformat.write(nb,ROOT/"Assignment5.ipynb")


if __name__=="__main__":
    s=Simulation();print(json.dumps(s.summary(),indent=2),flush=True)
    validation=run_validation(s)
    if not validation["all_passed"]:raise RuntimeError(validation)
    (OUT/"validation.json").write_text(json.dumps(validation,indent=2))
    (OUT/"summary.json").write_text(json.dumps(s.summary(),indent=2))
    export_tables(s); print("Tables and validation saved",flush=True)
    plot_results(s); print("Figures saved",flush=True)
    export_dashboard(s,validation);print("Offline interactive demo saved",flush=True)
    create_notebook();print("Report notebook saved",flush=True)
