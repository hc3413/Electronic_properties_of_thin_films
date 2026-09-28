'''Unit checks for the Hall analysis: carrier density, mobility and the sheet identity sigma = e*n*mu.

film_thickness = 1 is the repo's sentinel for a 2D sheet (R_H in m^2/C, n in cm^-2).
Any other film_thickness is a real thickness in m (R_H in m^3/C, n in cm^-3).

Run from the repo root:  MPLBACKEND=Agg EPOTFvenv/bin/python -m pytest tests
'''
import sys, io, contextlib
from pathlib import Path

import numpy as np
import pytest
import scipy.constants
from scipy.stats import linregress

# The analysis modules import each other by name from the repo root
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from Class_Import import PPMSData, import_ppms_data
from Functions_General import carrier_density, extract_ctf
from Functions_VDP import vdp_resistivity, vdp_hall
from Functions_HallBar import hallbar_resistivity, hallbar_hall

e = scipy.constants.e

# Raw DR012 run4 file used for the real-data check (skipped if the data folder is not available)
DR012_RUN4 = Path('/Users/horatiocox/Desktop/RUG_Postdoc/Experiments/Electrical/VDP_DR012_LSO_BSO_BSO_STO/Data/DR012_LaScO3_BaSnO3_BaSnO3__STO_silver_run4')


def mock_ppms(film_thickness, R_H, R_s, measurement_type='VDP'):
    '''PPMSData with a linear Hall resistance R_xy = (R_H / t) * B and a uniform sheet resistance R_s.
    For a 2D sheet t = film_thickness = 1, so R_xy = R_H * B.'''
    temps = np.array([10.0, 300.0])
    fields = np.linspace(-9, 9, 19)  # T, zero field in the middle as extract_ctf leaves it
    currents = np.linspace(-1e-5, 1e-5, 5)
    n_idx = 6 if measurement_type == 'VDP' else 4
    data = np.zeros((len(temps), len(fields), len(currents), 5, n_idx))
    data[..., 0, :] = temps[:, None, None, None]
    data[..., 1, :] = fields[None, :, None, None]
    data[..., 2, :] = currents[None, None, :, None]
    I = currents[None, None, :]
    R_xy = (R_H / film_thickness) * fields[None, :, None]
    if measurement_type == 'VDP':
        # Hall configurations are indices 0/1, symmetric square VDP configurations are 2-5 (R_A = R_B = R_s ln2 / pi)
        data[..., 4, 0] = R_xy * I
        data[..., 4, 1] = R_xy * I
        data[..., 4, 2:6] = (R_s * np.log(2) / np.pi * I)[..., None]
        hb_dimensions = None
    else:
        # Longitudinal arms are indices 0/1 (R = R_s * l / w), Hall arms are 2/3
        hb_dimensions = (30, 100)
        data[..., 4, 0] = R_s * hb_dimensions[1] / hb_dimensions[0] * I
        data[..., 4, 1] = R_s * hb_dimensions[1] / hb_dimensions[0] * I
        data[..., 4, 2] = R_xy * I
        data[..., 4, 3] = R_xy * I
    tf_av = np.stack(np.broadcast_arrays(temps[:, None], fields[None, :]), axis=-1)
    ppms = PPMSData(film_thickness=film_thickness, measurement_type=measurement_type, hb_dimensions=hb_dimensions,
                    carrier_type='n_type', filename='mock', sample_code='MOCK')
    ppms.data_np_nd = data
    ppms.ctf = [currents, temps, fields, len(currents), len(temps), len(fields)]
    ppms.tf_av = tf_av
    with contextlib.redirect_stdout(io.StringIO()):
        if measurement_type == 'VDP':
            vdp_hall(vdp_resistivity([ppms]))
        else:
            hallbar_hall(hallbar_resistivity([ppms]))
    return ppms


def sheet_case():
    '''2D sheet: n_s = 1e13 cm^-2, mu = 100 cm^2/Vs -> R_H = 62.42 ohm/T, R_s = 6242 ohm/sq'''
    n_s = 1.00e13
    mu = 100.0
    R_H = 1 / (e * n_s * 1e4)
    R_s = 1 / (e * n_s * 1e4 * mu * 1e-4)
    return n_s, mu, R_H, R_s


def film_case():
    '''3D film, t = 20 nm: n = 1e19 cm^-3, mu = 10 cm^2/Vs -> R_H = 6.242e-7 m^3/C, R_s = rho / t'''
    n, mu, t = 1.00e19, 10.0, 20e-9
    R_H = 1 / (e * n * 1e6)
    rho = 1 / (e * n * 1e6 * mu * 1e-4)
    return n, mu, t, R_H, rho / t


def test_carrier_density_helper():
    n_s, _, R_H, _ = sheet_case()
    assert carrier_density(R_H, 0.0, 1)[0] == pytest.approx(n_s, rel=1e-12)
    n, _, t, R_H, _ = film_case()
    assert carrier_density(R_H, 0.0, t)[0] == pytest.approx(n, rel=1e-12)
    # error propagation dn = dR_H / (R_H^2 e), with the same unit conversion as n
    assert carrier_density(R_H, 0.01 * R_H, t)[1] == pytest.approx(0.01 * n, rel=1e-12)


@pytest.mark.parametrize('measurement_type', ['VDP', 'HallBar'])
def test_round_trip_2d_sheet(measurement_type):
    n_s, mu, R_H, R_s = sheet_case()
    ppms = mock_ppms(1, R_H, R_s, measurement_type)
    hc = ppms.hall_coefficient[1]
    assert hc[5] == pytest.approx(R_H, rel=1e-9)   # R_H in ohm/T = m^2/C
    assert hc[7] == pytest.approx(n_s, rel=1e-9)   # n_s in cm^-2
    assert hc[9] == pytest.approx(mu, rel=1e-9)    # mu in cm^2/Vs
    assert 'cm^-2' in ppms.hall_coefficient_df.columns[7]
    assert ppms.hall_data_df.columns[2].startswith('R_xy')


@pytest.mark.parametrize('measurement_type', ['VDP', 'HallBar'])
def test_round_trip_3d_film(measurement_type):
    n, mu, t, R_H, R_s = film_case()
    ppms = mock_ppms(t, R_H, R_s, measurement_type)
    hc = ppms.hall_coefficient[1]
    assert hc[5] == pytest.approx(R_H, rel=1e-9)   # R_H in m^3/C
    assert hc[7] == pytest.approx(n, rel=1e-9)     # n in cm^-3
    assert hc[9] == pytest.approx(mu, rel=1e-9)    # mu in cm^2/Vs
    assert 'cm^-3' in ppms.hall_coefficient_df.columns[7]
    assert ppms.hall_data_df.columns[2].startswith('rho_xy')


@pytest.mark.parametrize('measurement_type', ['VDP', 'HallBar'])
def test_identity_sigma_e_n_mu(measurement_type):
    '''sigma = e*n*mu: e * n_s[cm^-2] * mu[cm^2/Vs] * R_s[ohm/sq] = 1 for a sheet,
    e * n[cm^-3] * mu[cm^2/Vs] * rho[ohm.cm] = 1 for a film'''
    n_s, mu, R_H, R_s = sheet_case()
    ppms = mock_ppms(1, R_H, R_s, measurement_type)
    Bi = int(ppms.ctf[5] / 2)
    X = e * ppms.hall_coefficient[:, 7] * ppms.hall_coefficient[:, 9] * ppms.res_data[:, Bi, 4]
    assert X == pytest.approx(1.0, rel=1e-9)

    n, mu, t, R_H, R_s = film_case()
    ppms = mock_ppms(t, R_H, R_s, measurement_type)
    rho_ohm_cm = ppms.res_data[:, Bi, 4] * 1e2
    X = e * ppms.hall_coefficient[:, 7] * ppms.hall_coefficient[:, 9] * rho_ohm_cm
    assert X == pytest.approx(1.0, rel=1e-9)


@pytest.mark.skipif(not DR012_RUN4.exists(), reason='DR012 raw data not available')
def test_dr012_real_data():
    '''DR012 run4 (Analysis_Publication settings): the identity holds at every temperature, and the
    302 K sheet density equals 1/(e R_H) from a direct fit of the raw R_xy vs B, in cm^-2.'''
    with contextlib.redirect_stdout(io.StringIO()):
        files = [f for f in import_ppms_data(DR012_RUN4.parent, film_thickness=1) if f.filename == DR012_RUN4.name]
        dat, _ = extract_ctf(files, ohm_m=1, single_zero_field=True)
        dat = vdp_hall(vdp_resistivity(dat), threshold=2)
    d = dat[0]
    Bi = int(d.ctf[5] / 2)
    X = e * d.hall_coefficient[:, 7] * d.hall_coefficient[:, 9] * d.res_data[:, Bi, 4]
    assert X == pytest.approx(1.0, rel=1e-9)

    # Direct fit of the raw file near 302 K: R_xy per field setpoint, then slope against B in tesla
    # Columns: T (K), field (Oe), source I (A), sense V (V), configuration index. 'New Temperature' lines are skipped as comments
    raw = np.genfromtxt(DR012_RUN4, skip_header=7, usecols=(0, 1, 2, 4, 5), comments='N')
    raw = raw[~np.isnan(raw).any(axis=1) & (raw[:, 0] > 300) & (raw[:, 0] < 304)]
    blocks = np.split(raw, np.flatnonzero(np.abs(np.diff(raw[:, 1])) > 50) + 1)
    B, R_xy = [], []
    for blk in blocks:
        R_A = linregress(blk[blk[:, 4] == 0, 2], blk[blk[:, 4] == 0, 3]).slope
        R_B = linregress(blk[blk[:, 4] == 1, 2], blk[blk[:, 4] == 1, 3]).slope
        B.append(blk[:, 1].mean() * 1e-4)  # Oe -> T
        R_xy.append((R_A + R_B) / 2)
    n_s_direct = 1e-4 / (e * linregress(B, R_xy).slope)  # m^-2 -> cm^-2
    Ti = int(np.argmin(np.abs(d.hall_coefficient[:, 0] - 302)))
    assert d.hall_coefficient[Ti, 7] == pytest.approx(n_s_direct, rel=1e-3)
    assert d.hall_coefficient[Ti, 7] == pytest.approx(4.86e14, rel=0.01)
