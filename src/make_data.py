"""Build real solar system CSV files from published orbital data.

Reads J2000 element tables embedded below plus an MPCORB file prefix
and writes heliocentric cartesian state vectors in AU, AU/day and
solar masses with header x,y,z,vx,vy,vz,m.
"""
import argparse
import csv
import math
import os
import sys

GAUSS_K = 0.01720209895
J2000_JD = 2451545.0
MSUN_KG = 1.98847e30
AU_KM = 149597870.7
DEG = math.pi / 180.0

PLANETS = [
    ("Mercury", 0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593, 1.660137e-7),
    ("Venus", 0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255, 2.447838e-6),
    ("Earth", 1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0, 3.0034896e-6),
    ("Mars", 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891, 3.227136e-7),
    ("Jupiter", 5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909, 9.5479187e-4),
    ("Saturn", 9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448, 2.8588564e-4),
    ("Uranus", 19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.96427630, 74.01692503, 4.366244e-5),
    ("Neptune", 30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574, 5.151389e-5),
]

DWARFS = [
    ("Pluto", 2457588.5, 39.6, 0.252, 17.1, 110.0, 114.0, 38.7, 0.00396, 6.553e-9),
    ("Ceres", 2461200.5, 2.77, 0.0797, 10.6, 80.2, 73.3, 274.0, 0.214, 4.719e-10),
    ("Eris", 2461200.5, 67.9, 0.438, 43.9, 36.0, 151.0, 212.0, 0.00176, 8.35e-9),
    ("Makemake", 2461200.5, 45.6, 0.159, 29.0, 79.3, 297.0, 170.0, 0.0032, 1.56e-9),
    ("Haumea", 2461200.5, 43.1, 0.194, 28.2, 122.0, 241.0, 223.0, 0.00349, 2.015e-9),
]

MOONS = [
    ("Moon", "Earth", 384400.0, 27.321661, 0.0549, 5.14, 0),
    ("Io", "Jupiter", 421800.0, 1.769138, 0.0, None, 0),
    ("Europa", "Jupiter", 671100.0, 3.551181, 0.0, None, 0),
    ("Ganymede", "Jupiter", 1070400.0, 7.154553, 0.0, None, 0),
    ("Callisto", "Jupiter", 1882700.0, 16.689018, 0.0, None, 0),
    ("Amalthea", "Jupiter", 181400.0, 0.498179, 0.0, None, 0),
    ("Mimas", "Saturn", 185500.0, 0.942422, 0.0, None, 0),
    ("Enceladus", "Saturn", 238000.0, 1.370218, 0.0, None, 0),
    ("Tethys", "Saturn", 294700.0, 1.887802, 0.0, None, 0),
    ("Dione", "Saturn", 377400.0, 2.736915, 0.0, None, 0),
    ("Rhea", "Saturn", 527100.0, 4.518212, 0.0, None, 0),
    ("Titan", "Saturn", 1221900.0, 15.945, 0.0, None, 0),
    ("Hyperion", "Saturn", 1481000.0, 21.276, 0.0, None, 0),
    ("Iapetus", "Saturn", 3560800.0, 79.3215, 0.0, None, 0),
    ("Miranda", "Uranus", 129900.0, 1.413479, 0.0, None, 0),
    ("Ariel", "Uranus", 191000.0, 2.520379, 0.0, None, 0),
    ("Umbriel", "Uranus", 266000.0, 4.144177, 0.0, None, 0),
    ("Titania", "Uranus", 436000.0, 8.705872, 0.0, None, 0),
    ("Oberon", "Uranus", 583500.0, 13.463234, 0.0, None, 0),
    ("Triton", "Neptune", 354800.0, 5.876854, 0.0, None, 1),
    ("Proteus", "Neptune", 117600.0, 1.122315, 0.0, None, 0),
    ("Nereid", "Neptune", 5513400.0, 360.136, 0.7507, None, 0),
    ("Charon", "Pluto", 19590.0, 6.387230, 0.0, None, 0),
    ("Nix", "Pluto", 48700.0, 24.8546, 0.0, None, 0),
    ("Hydra", "Pluto", 64700.0, 38.2018, 0.0, None, 0),
]

MOON_MASS_KG = {
    "Moon": 7.342e22,
    "Io": 8.932e22,
    "Europa": 4.8e22,
    "Ganymede": 1.482e23,
    "Callisto": 1.076e23,
    "Amalthea": 2.08e18,
    "Mimas": 3.75e19,
    "Enceladus": 1.08e20,
    "Tethys": 6.17e20,
    "Dione": 1.096e21,
    "Rhea": 2.31e21,
    "Titan": 1.345e23,
    "Hyperion": 5.62e18,
    "Iapetus": 1.81e21,
    "Miranda": 6.59e19,
    "Ariel": 1.35e21,
    "Umbriel": 1.17e21,
    "Titania": 3.53e21,
    "Oberon": 3.01e21,
    "Triton": 2.14e22,
    "Proteus": 5.0e19,
    "Nereid": 3.1e19,
    "Charon": 1.586e21,
    "Nix": 4.5e16,
    "Hydra": 4.8e16,
}


def kepler_E(M, e):
    E = M
    for _ in range(50):
        delta = (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
        E -= delta
        if abs(delta) < 1e-14:
            break
    return E


def elements_to_cartesian(a, e, inc, Om, om, M):
    E = kepler_E(M, e)
    cE = math.cos(E)
    sE = math.sin(E)
    sq = math.sqrt(1.0 - e * e)
    xp = a * (cE - e)
    yp = a * sq * sE
    n = GAUSS_K / (a * math.sqrt(a))
    denom = 1.0 - e * cE
    vxp = -a * n * sE / denom
    vyp = a * sq * n * cE / denom
    cO = math.cos(Om)
    sO = math.sin(Om)
    co = math.cos(om)
    so = math.sin(om)
    ci = math.cos(inc)
    si = math.sin(inc)
    x = (cO * co - sO * so * ci) * xp + (-cO * so - sO * co * ci) * yp
    y = (sO * co + cO * so * ci) * xp + (-sO * so + cO * co * ci) * yp
    z = (so * si) * xp + (co * si) * yp
    vx = (cO * co - sO * so * ci) * vxp + (-cO * so - sO * co * ci) * vyp
    vy = (sO * co + cO * so * ci) * vxp + (-sO * so + cO * co * ci) * vyp
    vz = (so * si) * vxp + (co * si) * vyp
    return (x, y, z, vx, vy, vz)


def state_from_elements(name, a, e, inc, node, argp, M, mass):
    pos = elements_to_cartesian(a, e, math.radians(inc), math.radians(node), math.radians(argp), math.radians(M))
    return name, pos, mass, math.radians(inc), math.radians(node)


def planet_state(entry):
    name, a, e, inc, L, varpi, node, mass = entry
    return state_from_elements(name, a, e, inc, node, varpi - node, (L - varpi) % 360.0, mass)


def dwarf_state(entry):
    name, epoch, a, e, inc, node, peri, M, n, mass = entry
    return state_from_elements(name, a, e, inc, node, peri, (M + n * (J2000_JD - epoch)) % 360.0, mass)


def moon_state(entry, parent_pos, parent_inc, parent_node, phase):
    name, _parent, a_km, period, e, tilt, retro = entry
    a = a_km / AU_KM
    E = kepler_E(phase, e)
    cE = math.cos(E)
    sE = math.sin(E)
    sq = math.sqrt(1.0 - e * e)
    xp = a * (cE - e)
    yp = a * sq * sE
    rate = 2.0 * math.pi / period
    denom = 1.0 - e * cE
    vxp = -a * rate * sE / denom
    vyp = a * sq * rate * cE / denom
    if retro:
        yp = -yp
        vyp = -vyp
    if tilt is None:
        inc = parent_inc
        node = parent_node
    else:
        inc = math.radians(tilt)
        node = parent_node
    ci = math.cos(inc)
    si = math.sin(inc)
    cO = math.cos(node)
    sO = math.sin(node)
    y1 = yp * ci
    z1 = yp * si
    vy1 = vyp * ci
    vz1 = vyp * si
    x = xp * cO - y1 * sO
    y = xp * sO + y1 * cO
    z = z1
    vx = vxp * cO - vy1 * sO
    vy = vxp * sO + vy1 * cO
    vz = vz1
    mass = MOON_MASS_KG[name] / MSUN_KG
    return (parent_pos[0] + x, parent_pos[1] + y, parent_pos[2] + z,
            parent_pos[3] + vx, parent_pos[4] + vy, parent_pos[5] + vz, mass)


def build_core_rows():
    states = {}
    planes = {}
    rows = [(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]
    for entry in PLANETS:
        name, pos, mass, inc, node = planet_state(entry)
        states[name] = pos
        planes[name] = (inc, node)
        rows.append((pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], mass))
    for entry in DWARFS:
        name, pos, mass, inc, node = dwarf_state(entry)
        states[name] = pos
        planes[name] = (inc, node)
        rows.append((pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], mass))
    siblings = {}
    for entry in MOONS:
        siblings.setdefault(entry[1], []).append(entry)
    for parent, group in siblings.items():
        for index, entry in enumerate(group):
            phase = 2.0 * math.pi * index / len(group)
            inc, node = planes[parent]
            rows.append(moon_state(entry, states[parent], inc, node, phase))
    return rows


def mass_from_H(H):
    diameter_km = 1329.0 / math.sqrt(0.15) * 10.0 ** (-H / 5.0)
    radius_m = diameter_km * 500.0
    volume = 4.0 / 3.0 * math.pi * radius_m ** 3
    return volume * 2500.0 / MSUN_KG


def parse_mpcorb(path, limit):
    out = []
    with open(path, errors="ignore") as f:
        for line in f:
            if len(line) < 100:
                continue
            try:
                a = float(line[92:103].strip())
                e = float(line[70:79].strip())
                inc = float(line[59:68].strip())
                Om = float(line[48:57].strip())
                om = float(line[37:46].strip())
                M = float(line[26:35].strip())
                H = float(line[8:13].strip())
            except ValueError:
                continue
            if not (0.0 <= e < 1.0):
                continue
            if not (1.8 < a < 3.6):
                continue
            out.append((a, e, math.radians(inc), math.radians(Om), math.radians(om), math.radians(M), H))
            if len(out) >= limit:
                break
    return out


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["x", "y", "z", "vx", "vy", "vz", "m"])
        for row in rows:
            writer.writerow(["%.12e" % v for v in row])


def main(argv):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-out", default=os.path.join(root, "data", "solar_system.csv"))
    parser.add_argument("--mpcorb", default=None)
    parser.add_argument("--limit", type=int, default=8000)
    parser.add_argument("--asteroids-out", default=os.path.join(root, "data", "asteroids.csv"))
    args = parser.parse_args(argv)
    core_rows = build_core_rows()
    write_csv(args.core_out, core_rows)
    print("core bodies: " + str(len(core_rows)) + " -> " + args.core_out)
    if args.mpcorb:
        elements = parse_mpcorb(args.mpcorb, args.limit)
        rows = []
        for a, e, inc, Om, om, M, H in elements:
            x, y, z, vx, vy, vz = elements_to_cartesian(a, e, inc, Om, om, M)
            rows.append((x, y, z, vx, vy, vz, mass_from_H(H)))
        write_csv(args.asteroids_out, rows)
        print("asteroids: " + str(len(rows)) + " -> " + args.asteroids_out)


if __name__ == "__main__":
    main(sys.argv[1:])
