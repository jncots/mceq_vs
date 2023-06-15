from pathlib import Path
import h5py
import numpy as np
from particle import Particle


all_particles_dict = {p.pdgid: p for p in Particle.findall()}
xdepth_list = np.array([21, 86, 210, 492, 1033], dtype=np.float32)


def corsika_hist_en(en_bins, h5file):
    energy_hist = dict()
    print(f"Read file:\n.../{h5file.relative_to(h5file.parents[2])}")
    with h5py.File(h5file, "r") as corsika_data:
        
        num_primaries = corsika_data["num_primaries"]
        print(f"Number of primaries = {np.array(num_primaries)*1.0:e}")
        
        pdg_dict = {}
        for pdg in corsika_data:
            try:
                pdg_dict[int(pdg)] = [int(ixdepth) for ixdepth in corsika_data[pdg].keys()]
            except:
                pass
        
        for pdg in pdg_dict:
            xd_dict = []
            for ixdepth in pdg_dict[pdg]:
                mcdata = np.array(corsika_data[str(pdg)][str(ixdepth)]["energy [GeV]"])
                hist, bin_edges = np.histogram(mcdata, bins = en_bins)
                num_tot = np.sum(hist)
                xd_dict.append((hist, bin_edges, xdepth_list[ixdepth], np.array(num_primaries)*1.0, num_tot))        
            energy_hist[pdg] = (xd_dict, f"${all_particles_dict[pdg].latex_name}$")
    return energy_hist


def combined_data_en(energy_hist, pdgs, ixdepth):
    
    hist = None
    label = ""
    for pdg in pdgs:
        
        if hist is None:
            hist = energy_hist[pdg][0][ixdepth][0]
        else:
            hist = energy_hist[pdg][0][ixdepth][0] + hist    
        
        label += f"+{energy_hist[pdg][1]}"
    
    label_p = f"{label[1:]}"
    label_x = f"X={energy_hist[pdg][0][ixdepth][2]}"
    
    return (hist, energy_hist[pdg][0][ixdepth][1], 
            label_p, label_x, 
            energy_hist[pdg][0][ixdepth][3], # num_primaries
            )


def corsika_en_theta_2dhist(en_bins, theta_bins, h5file):
    
    base_dir = Path("/hetghome/antonpr/xmax_sigma/data_comparison")
    h5file = base_dir / h5file
    
    energy_hist = dict()
    with h5py.File(h5file, "r") as corsika_data:
        num_primaries = corsika_data["num_primaries"]
        
        pdg_dict = {}
        for pdg in corsika_data:
            try:
                pdg_dict[int(pdg)] = [int(ixdepth) for ixdepth in corsika_data[pdg].keys()]
            except:
                pass
        
        for pdg in pdg_dict:
            print(f"pdg = {pdg}")
            ix = [int(key) for key in corsika_data[str(pdg)].keys()]
            xd_dict = []
            for i, xdepth in enumerate(xdepth_list):
                mc_energy = np.array(corsika_data[str(pdg)][str(ix[i])]["energy [GeV]"])
                mc_theta = np.array(corsika_data[str(pdg)][str(ix[i])]["theta [rad]"])
                print(f"xdepth={xdepth}, number={len(mc_energy)*1e0:.4e}")
                
                # If Omega is used for binning
                # mc_omega = 2*np.pi*(1 - np.cos(mc_theta))
                # mc_theta = mc_omega
                
                hist, en_edges, theta_edges = np.histogram2d(mc_energy, mc_theta, 
                                                             bins = [en_bins, theta_bins])
                
                # If Omega is used for binning
                # omega_edges = np.arccos(1 - theta_edges/(2*np.pi))
                # theta_edges = omega_edges
                
                # if pdg in [-13, 13]:
                #     print(f"Hist = {hist}")
                #     print(f"Hist = {np.sqrt(hist)/hist}")
                #     print(f"Xdepth = {xdepth}, pdg = ${pdg_dict[pdg].latex_name}$")
                #     # print(f"Hist = {theta_edges}")
                
                xd_dict.append((hist/num_primaries, en_edges, theta_edges, xdepth, np.sqrt(hist)/num_primaries))         
            energy_hist[pdg] = (xd_dict, f"${pdg_dict[pdg].latex_name}$")
    return energy_hist


def combined_ang_data(energy_hist, pdgs):    
    dist_xdepth = []
    for ixdepth in range(3):
        
        dist_en = []
        for ind_energy in range(energy_hist[13][0][ixdepth][1].size - 1):
            
            hist = None
            hist_error = None
            label = ""
            for pdg in pdgs:    
                ang_dist = energy_hist[pdg][0][ixdepth][0][ind_energy]
                ang_dist_error = energy_hist[pdg][0][ixdepth][4][ind_energy]
                en_bins = energy_hist[pdg][0][ixdepth][1]
                cur_en_label = f"[{en_bins[ind_energy]}, {en_bins[ind_energy + 1]}]"
                cur_en = np.array([en_bins[ind_energy], en_bins[ind_energy + 1]])
                ang_bins = energy_hist[pdg][0][ixdepth][2]
                xdepth = energy_hist[pdg][0][ixdepth][3]
        
                if hist is None:
                    hist = ang_dist
                else:
                    hist = ang_dist + hist   

                if hist_error is None:
                    hist_error = ang_dist_error**2
                else:
                    hist_error = ang_dist_error**2 + hist_error   
                
                # print(f"pdg = {pdg}, sum = {np.sum(hist)}, size = {hist.size}")


                
                label += f"+{energy_hist[pdg][1]}"
    
                label_p = f"{label[1:]}"
                label_x = f"{xdepth}"

            hist_error = np.sqrt(hist_error)
            # print(f"hist_error = {hist_error/hist}")
            # print(f"ang_dist_error = {ang_dist_error/ang_dist}")

            dist_en.append((hist, ang_bins, cur_en_label, label_x, label_p, cur_en, hist_error))
        dist_xdepth.append(dist_en)    
    return dist_xdepth

def print_h5file_structure(h5file):
    """Print a structure of db
    Args:
        h5file (pathlib.Path): path to *.h5 file
    """
    import nexusformat.nexus as nx
    f = nx.nxload(h5file)
    print(f.tree)
    
    # with h5py.File(h5file, "r") as corsika_data:
    #     num_primaries = corsika_data["num_primaries"]
    #     print(np.array(num_primaries))
      
    
if __name__ == "__main__":
    import nexusformat.nexus as nx
    base_dir = Path("/hetghome/antonpr/projects/mceq_vs/data/corsika")
    h5file = base_dir/"03_single_muon_shower/01_single_muon_run/01_single_muon.h5"
    print_h5file_structure(h5file)
    
        