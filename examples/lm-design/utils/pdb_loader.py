# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# Portions of this file were adapted from the open source code
# https://github.com/RosettaCommons/trRosetta2/tree/main
# Modifications were made to loader() and get_coords6d() functions:
# 1. addition of `set_diagonal` flag.
# 2. addition of `allow_missing_residue_coords` flag.
# Original License information below.

# MIT License

# Copyright (c) 2021 Ivan Anishchenko, Minkyung Baek, Naozumi Hiranuma, Ian Humphrey

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import numpy as np
import scipy
import scipy.spatial
import logging

logger = logging.getLogger(__name__)


PDB_LOADER_PARAMS_DEFAULT = {  # trRosetta-v2
    "DMIN"    : 2.0,
    "DMAX"    : 20.0,
    "DBINS"   : 36,
    "ABINS"   : 36,
    "WMIN"    : 0.8,
    "LMIN"    : 150,
    "LMAX"    : 400,
    "EPOCHS"  : 10,
    "NCPU"    : 8,
    "SLICE"   : "CONT",
    "contact_bin_cutoff": (0, 11) # both inclusive
}


# calculate dihedral angles defined by 4 sets of points
def get_dihedrals(a, b, c, d):

    b0 = -1.0*(b - a)
    b1 = c - b
    b2 = d - c

    b1 /= np.linalg.norm(b1, axis=-1)[:,None]

    v = b0 - np.sum(b0*b1, axis=-1)[:,None]*b1
    w = b2 - np.sum(b2*b1, axis=-1)[:,None]*b1

    x = np.sum(v*w, axis=-1)
    y = np.sum(np.cross(b1, v)*w, axis=-1)

    return np.arctan2(y, x)


# calculate planar angles defined by 3 sets of points
def get_angles(a, b, c):

    v = a - b
    v /= np.linalg.norm(v, axis=-1)[:,None]

    w = c - b
    w /= np.linalg.norm(w, axis=-1)[:,None]

    x = np.sum(v*w, axis=1)

    return np.arccos(x)


def nonans(arr):
    return not (arr != arr).any()


# get 6d coordinates from x,y,z coords of N,Ca,C atoms
def get_coords6d(xyz, dmax, allow_missing_residue_coords=False):

    nres = xyz.shape[1]
    # three anchor atoms
    N  = xyz[0]
    Ca = xyz[1]
    C  = xyz[2]
    if not nonans(xyz):
        assert allow_missing_residue_coords, "Missing residue coordinates"
        logger.warning("PDB xyz contains NaNs!! Loading anyway...")
    # recreate Cb given N,Ca,C
    b = Ca - N
    c = C - Ca
    a = np.cross(b, c)
    Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca

    # fast neighbors search to collect all
    # Cb-Cb pairs within dmax
    kdCb = scipy.spatial.cKDTree(Cb)
    indices = kdCb.query_ball_tree(kdCb, dmax)

    # indices of contacting residues
    idx = np.array([[i,j] for i in range(len(indices)) for j in indices[i] if i != j]).T
    idx0 = idx[0]
    idx1 = idx[1]

    # Cb-Cb distance matrix
    dist6d = np.full((nres, nres), 999.9)
    dist6d[idx0,idx1] = np.linalg.norm(Cb[idx1]-Cb[idx0], axis=-1)

    # matrix of Ca-Cb-Cb-Ca dihedrals
    omega6d = np.zeros((nres, nres))
    omega6d[idx0,idx1] = get_dihedrals(Ca[idx0], Cb[idx0], Cb[idx1], Ca[idx1])

    # matrix of polar coord theta
    theta6d = np.zeros((nres, nres))
    theta6d[idx0,idx1] = get_dihedrals(N[idx0], Ca[idx0], Cb[idx0], Cb[idx1])

    # matrix of polar coord phi
    phi6d = np.zeros((nres, nres))
    phi6d[idx0,idx1] = get_angles(Ca[idx0], Cb[idx0], Cb[idx1])

    if nonans(xyz):
        # Set dist values to nan where there's no pdb 3d coordinates data
        nan_residues = ((Cb != Cb).sum(-1) > 0)
        nan_residues_mask = (nan_residues[None] | nan_residues[None].T)
        dist6d[nan_residues_mask] = np.nan
        omega6d[nan_residues_mask] = np.nan
        theta6d[nan_residues_mask] = np.nan
        phi6d[nan_residues_mask] = np.nan

    return dist6d, omega6d, theta6d, phi6d


def parse_PDB(pdb_path, atoms=['N','CA','C'], chain=None, return_aligned_seq=False):
    '''
    input:  pdb_path = PDB filename
            atoms = atoms to extract (optional)
    output: (length, atoms, coords=(x,y,z)), sequence
    '''
    xyz,seq,doubles,min_resn,max_resn = {},{},{},np.inf,-np.inf
    with open(pdb_path,"rb") as f:
        for line in f:
            line = line.decode("utf-8","ignore").rstrip()

            if line[:6] == "HETATM" and line[17:17+3] == "MSE":
                line = line.replace("HETATM","ATOM  ")
                line = line.replace("MSE","MET")

            if line[:4] == "ATOM":
                ch = line[21:22]
                if ch == chain or chain is None:
                    atom = line[12:12+4].strip()
                    resi = line[17:17+3]
                    resi_extended = line[16:17 + 3].strip()
                    resn = line[22:22+5].strip()
                    x,y,z = [float(line[i:(i+8)]) for i in [30,38,46]]

                    if resn[-1].isalpha(): resa,resn = resn[-1],int(resn[:-1])-1
                    else: resa,resn = "",int(resn)-1
                    if resn < min_resn: min_resn = resn
                    if resn > max_resn: max_resn = resn
                    if resn not in xyz: xyz[resn] = {}
                    if resa not in xyz[resn]: xyz[resn][resa] = {}
                    if resn not in seq: seq[resn] = {}
                    if resa not in seq[resn]:
                        seq[resn][resa] = resi
                    elif seq[resn][resa] != resi_extended:
                        # doubles mark locations in the pdb file where multi residue entries are
                        # present. There's a known bug in TmAlign binary that doesn't read / skip
                        # these entries, so we mark them to create a sequence that is aligned with
                        # gap tokens in such locations.
                        doubles[resn] = True

                    if atom not in xyz[resn][resa]:
                        xyz[resn][resa][atom] = np.array([x,y,z])

    # convert to numpy arrays, fill in missing values
    seq_,xyz_,aligned_seq_ = [],[],[]
    for resn in range(min_resn,max_resn+1):
        if resn in seq:
            for k in sorted(seq[resn]):
                seq_.append(aa_3_N.get(seq[resn][k],20))
                aligned_seq_.append(seq_[-1]) if resn not in doubles else aligned_seq_.append(20)
        else:
            seq_.append(20)
        if resn in xyz:
            for k in sorted(xyz[resn]):
                for atom in atoms:
                    if atom in xyz[resn][k]: xyz_.append(xyz[resn][k][atom])
                    else: xyz_.append(np.full(3,np.nan))
        else:
            for atom in atoms: xyz_.append(np.full(3,np.nan))

    res = [np.array(xyz_).reshape(-1,len(atoms),3), N_to_AA(np.array(seq_))]
    if return_aligned_seq:
        res += [N_to_AA(np.array(aligned_seq_))]
    return res

def N_to_AA(x):
    # [[0,1,2,3]] -> ["ARND"]
    x = np.array(x);
    if x.ndim == 1: x = x[None]
    return ["".join([aa_N_1.get(a,"-") for a in y]) for y in x]

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','GAP']

aa_3_N = {a:n for n,a in enumerate(alpha_3)}
aa_N_1 = {n:a for n,a in enumerate(alpha_1)}

def loader(pdb_path, params, chain=None, set_diagonal=False, allow_missing_residue_coords=False):
    """
    Args:
        pdb_path: path to pdb file
        params: dict with hyperparams
        chain: <unused>
        set_diagonal: sets diagonal to specific default values.
    """
    orig_xyz, seq = parse_PDB(pdb_path, atoms=['N','CA','C'], chain=None)
    xyz = np.transpose(orig_xyz, (1,0,2))
    idx = np.arange(xyz.shape[1])
    # get 6D coords
    d,o,t,p = get_coords6d(xyz, params['DMAX'], allow_missing_residue_coords)
    no_angles = any([np.isnan(x).any() for x in [o,t,p]])

    # bin 6D coords
    if params['DMIN'] == 2.5:
        dbins = np.linspace(params['DMIN'], params['DMAX'], params['DBINS'])
    else:
        dstep = (params['DMAX'] - params['DMIN']) / params['DBINS']
        dbins = np.linspace(params['DMIN'] + dstep, params['DMAX'], params['DBINS'])

    angle_step = 2.0 * np.pi / params['ABINS']
    ab360 = np.linspace(-np.pi + angle_step, np.pi, params['ABINS'])
    phi_last_bin = int(params['ABINS']/2)
    if "PHI_BINS" in params:
        phi_last_bin = params['PHI_BINS']
    ab180 = np.linspace(angle_step, np.pi, phi_last_bin)

    db = np.digitize(d, dbins).astype(np.uint8)  # distance
    ob = np.digitize(o, ab360).astype(np.uint8)  # omega
    tb = np.digitize(t, ab360).astype(np.uint8)  # theta
    pb = np.digitize(p, ab180).astype(np.uint8)  # phi

    # synchronize 'no contact' bins
    ob[db == params['DBINS']] = params['ABINS']
    tb[db == params['DBINS']] = params['ABINS']
    pb[db == params['DBINS']] = phi_last_bin

    if set_diagonal:
        db[np.eye(db.shape[0]).astype(bool)] = 0
        ob[np.eye(ob.shape[0]).astype(bool)] = int(params['ABINS']/2)
        tb[np.eye(tb.shape[0]).astype(bool)] = int(params['ABINS']/2)
        pb[np.eye(pb.shape[0]).astype(bool)] = 0

        # (11/20/2021)
        # Added masking of this diagonal as well.
        # This matches the way contact predictions were trained.
        d[np.eye(db.shape[0]).astype(bool)] = 0

    # stack all coords together
    if no_angles:
        c6d = db[None, :]  # only distogram, unsqueezed(0)
    else:
        c6d = np.stack([db,ob,tb,pb], axis=0)

    # slice long chains
    L = idx.shape[0]
    start,stop,nres = 0,L,L
    sel = np.arange(0, L)
    if L > params['LMAX']:

        if params['SLICE'] == 'CONT':
            # slice continuously
            nres = np.random.randint(params['LMIN'], params['LMAX'])
            logger.warning("Slicing long chain %s into %s residues", pdb_path, nres)
#            nres = params['LMAX']
            start = np.random.randint(L-nres+1)
#            start = 0
            stop = start + nres
            sel = np.arange(start,stop)

    if (idx[sel][1:]-idx[sel][:-1]==1).all(): #check if idx has no gaps
        feed_dict = {
            "idx"      : idx[sel],
            "idx_raw"  : idx,
            "coords6d" : c6d[:,sel,:][:,:,sel],
            "sel"      : sel,
            "dist"     : d[sel,:][:,sel],
            "fullseq"  : [seq[0]],
            'no_angles': no_angles,
            'xyz'      : orig_xyz
        }
    else:
        feed_dict = {
            "idx"      : idx[sel],
            "idx_raw"  : idx,
            "coords6d" : np.zeros(1),
            "sel"      : sel,
            "dist"     : d[sel,:][:,sel],
            "fullseq"  : [seq[0]],
            'no_angles': no_angles,
            'xyz'      : orig_xyz
        }
    return feed_dict
