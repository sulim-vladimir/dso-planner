#!/usr/bin/env python3
"""
dso_window.py - find NGC / IC / Messier objects that pass through a
"window" of sky (an azimuth range and altitude range) during a time period.

Default window: azimuth NW..N (315..360 deg), altitude 10..45 deg.

Example:
    ./dso_window.py --lat 50.45 --lon 30.52 --date 2026-11-01 --from 21:00 --to 01:00
    ./dso_window.py --date 2026-11-01 --from 21:00 --to 23:00 --types galaxy --max-mag 11

Catalog data: OpenNGC by Mattia Verga (CC-BY-SA-4.0), files in ./data/
"""

import argparse
import csv
import json
import math
import os
import re
import sys
from datetime import date, datetime, timedelta, timezone

try:
    from zoneinfo import ZoneInfo
except ImportError:  # Python < 3.9
    ZoneInfo = None

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
CONFIG_FILE = os.path.join(HERE, "config.json")
FAVORITES_FILE = os.path.join(HERE, "favorites.txt")


def norm_id(s):
    """'M 101', 'm101', 'NGC0224' -> 'M101', 'M101', 'NGC224' (for matching names)."""
    s = re.sub(r"\s+", "", s).upper()
    m = re.match(r"^(NGC|IC|M)0*(\d+.*)$", s)
    return m.group(1) + m.group(2) if m else s


def load_favorites(path):
    """Read favorites file: one object per line, '#' starts a comment."""
    if not path or not os.path.exists(path):
        return set()
    with open(path, encoding="utf-8") as f:
        return {norm_id(line.split("#")[0]) for line in f if line.split("#")[0].strip()}


def is_favorite(o, favs):
    return bool(favs) and any(norm_id(n) in favs for n in [o["name"], *o["aliases"], o["common"]] if n)

# ---------------------------------------------------------------------------
# Object types
# ---------------------------------------------------------------------------
# OpenNGC type code -> (group, readable name)
TYPES = {
    "G": ("galaxy", "Galaxy"),
    "GPair": ("galaxy", "Galaxy pair"),
    "GTrpl": ("galaxy", "Galaxy triplet"),
    "GGroup": ("galaxy", "Galaxy group"),
    "OCl": ("cluster", "Open cluster"),
    "GCl": ("cluster", "Globular cluster"),
    "*Ass": ("cluster", "Stellar association"),
    "Cl+N": ("cluster", "Cluster + nebula"),
    "PN": ("nebula", "Planetary nebula"),
    "HII": ("nebula", "HII region"),
    "EmN": ("nebula", "Emission nebula"),
    "RfN": ("nebula", "Reflection nebula"),
    "SNR": ("nebula", "Supernova remnant"),
    "Neb": ("nebula", "Nebula"),
    "DrkN": ("nebula", "Dark nebula"),
    "*": ("other", "Star"),
    "**": ("other", "Double star"),
    "Nova": ("other", "Nova"),
    "Other": ("other", "Other"),
}

COMPASS = {
    "N": 0, "NNE": 22.5, "NE": 45, "ENE": 67.5, "E": 90, "ESE": 112.5,
    "SE": 135, "SSE": 157.5, "S": 180, "SSW": 202.5, "SW": 225,
    "WSW": 247.5, "W": 270, "WNW": 292.5, "NW": 315, "NNW": 337.5,
}


def parse_az(text):
    """Accept a number of degrees or a compass point (N, NW, NNW ...)."""
    t = str(text).strip().upper()
    if t in COMPASS:
        return float(COMPASS[t])
    return float(t) % 360.0


def pretty_name(code):
    """NGC0224 -> NGC 224, IC0001 -> IC 1, M045 -> M45."""
    m = re.match(r"^(NGC|IC|M)0*(\d+)(.*)$", code)
    return f"{m.group(1)} {m.group(2)}{m.group(3)}" if m else code


def hms_to_deg(s):
    h, m, sec = s.split(":")
    return (float(h) + float(m) / 60 + float(sec) / 3600) * 15.0


def dms_to_deg(s):
    sign = -1.0 if s.strip().startswith("-") else 1.0
    d, m, sec = s.strip().lstrip("+-").split(":")
    return sign * (float(d) + float(m) / 60 + float(sec) / 3600)


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Catalog loading
# ---------------------------------------------------------------------------
def load_catalog():
    objects = []
    files = [("NGC.csv", False), ("addendum.csv", True)]
    for fname, addendum in files:
        path = os.path.join(DATA_DIR, fname)
        if not os.path.exists(path):
            if addendum:
                continue
            sys.exit(f"Catalog file not found: {path}\nRun ./update_catalog.sh to download it.")
        with open(path, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter=";"):
                t = r["Type"]
                if t in ("Dup", "NonEx") or not r["RA"] or not r["Dec"]:
                    continue
                m_num = int(r["M"]) if r["M"] else None
                # The addendum holds many non-NGC/IC objects; keep only Messier ones.
                if addendum and m_num is None:
                    continue
                name = r["Name"]
                if name.startswith("NGC"):
                    cat = "ngc"
                elif name.startswith("IC"):
                    cat = "ic"
                else:
                    cat = "other"

                vmag, bmag = fnum(r["V-Mag"]), fnum(r["B-Mag"])
                if vmag is not None:
                    mag, mag_band = vmag, "V"
                elif bmag is not None:
                    mag, mag_band = bmag, "B"
                else:
                    mag, mag_band = None, ""

                maj, mnr = fnum(r["MajAx"]), fnum(r["MinAx"])
                if maj is not None and mnr is None:
                    mnr = maj

                # Surface brightness (mag/arcsec^2): catalog value when present,
                # otherwise estimated from magnitude and size (marked with "~").
                sb, sb_est = fnum(r["SurfBr"]), False
                if sb is None and mag is not None and maj:
                    area = math.pi / 4 * (maj * 60) * (mnr * 60)
                    sb, sb_est = mag + 2.5 * math.log10(area), True

                group, tname = TYPES.get(t, ("other", t))
                if t == "G" and r["Hubble"]:
                    tname = f"Galaxy {r['Hubble']}"

                aliases = []
                if m_num:
                    aliases.append(f"M {m_num}")
                if cat != "ngc" and r["NGC"]:
                    aliases += [f"NGC {int(x)}" if x.isdigit() else f"NGC {x}"
                                for x in r["NGC"].split(",")]
                if cat != "ic" and r["IC"]:
                    aliases += [f"IC {int(x)}" if x.isdigit() else f"IC {x}"
                                for x in r["IC"].split(",")]

                objects.append({
                    "name": pretty_name(name),
                    "messier": m_num,
                    "cat": cat,
                    "aliases": aliases,
                    "common": r["Common names"].split(",")[0] if r["Common names"] else "",
                    "code": t,
                    "group": group,
                    "type": tname,
                    "mag": mag,
                    "mag_band": mag_band,
                    "sb": sb,
                    "sb_est": sb_est,
                    "maj": maj,
                    "min": mnr,
                    "pa": fnum(r["PosAng"]),
                    "ra_s": r["RA"],
                    "dec_s": r["Dec"],
                    "ra": hms_to_deg(r["RA"]),
                    "dec": dms_to_deg(r["Dec"]),
                    "const": r["Const"],
                })
    return objects


GREEK = {"Alp": "α", "Bet": "β", "Gam": "γ", "Del": "δ", "Eps": "ε", "Zet": "ζ", "Eta": "η", "The": "θ",
         "Iot": "ι", "Kap": "κ", "Lam": "λ", "Mu": "μ", "Nu": "ν", "Xi": "ξ", "Omi": "ο", "Pi": "π",
         "Rho": "ρ", "Sig": "σ", "Tau": "τ", "Ups": "υ", "Phi": "φ", "Chi": "χ", "Psi": "ψ", "Ome": "ω"}


def load_constellations():
    """Constellation stick figures (Stellarium, western sky culture) with star positions from HYG.
    Returns (stars, figures): stars = {hip: (ra, dec, mag, label)}, figures = [(name, [(hip, hip), ...])],
    or None when the data files are missing."""
    paths = [os.path.join(DATA_DIR, f) for f in
             ("constellation_stars.csv", "constellation_lines.fab", "constellation_names.fab")]
    if not all(os.path.exists(p) for p in paths):
        return None
    stars = {}
    with open(paths[0], newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            b = r["bayer"].split("-")
            greek = GREEK.get(b[0], b[0]) + (b[1] if len(b) > 1 else "")
            desig = f"{greek or r['flam']} {r['con']}".strip()
            label = " · ".join(x for x in (r["name"], desig) if x)
            stars[int(r["hip"])] = (float(r["ra"]), float(r["dec"]), float(r["mag"]), label)
    names = {}
    with open(paths[2], encoding="utf-8") as f:
        for ln in f:
            p = ln.split("\t")
            if len(p) > 1:
                names[p[0].strip()] = p[1].strip().strip('"')
    figures = []
    with open(paths[1], encoding="utf-8") as f:
        for ln in f:
            p = ln.split()
            if len(p) < 4:
                continue
            hips = [int(h) for h in p[2:]]
            segs = [(a, b) for a, b in zip(hips[0::2], hips[1::2]) if a in stars and b in stars]
            figures.append((names.get(p[0], p[0]), segs))
    return stars, figures


def constellation_payload(jd, lst_deg, lat, win_center):
    """Alt/az of the constellation figures at one moment, limited to the part of the sky
    around the window that the report's sky view can show."""
    data = load_constellations()
    if data is None:
        return None
    stars, figures = data
    hips = list(stars)
    ra, dec = precess(np.array([stars[h][0] for h in hips]), np.array([stars[h][1] for h in hips]), jd)
    alt, az = alt_az(ra[:, None], dec[:, None], lat, np.array([lst_deg]))
    pos = {h: (float(alt[i, 0]), float(az[i, 0])) for i, h in enumerate(hips)}

    def near(h):
        a, z = pos[h]
        return a > -10 and abs((z - win_center + 540) % 360 - 180) < 110

    out_stars, index, out_figs = [], {}, []
    for name, segs in figures:
        keep = [(a, b) for a, b in segs if near(a) or near(b)]
        if not keep:
            continue
        for h in {h for seg in keep for h in seg}:
            if h not in index:
                index[h] = len(out_stars)
                a, z = pos[h]
                out_stars.append([round(a, 2), round(z, 2), stars[h][2], stars[h][3]])
        out_figs.append({"name": name, "seg": [[index[a], index[b]] for a, b in keep]})
    return {"stars": out_stars, "figures": out_figs}


# ---------------------------------------------------------------------------
# Astronomy (accurate to a few arcminutes - plenty for choosing targets)
# ---------------------------------------------------------------------------
def julian_date(dt_utc):
    return dt_utc.timestamp() / 86400.0 + 2440587.5


def precess(ra_deg, dec_deg, jd):
    """Precess J2000 coordinates to the date (IAU 1976)."""
    T = (jd - 2451545.0) / 36525.0
    a2r = math.pi / (180 * 3600)
    zeta = (2306.2181 * T + 0.30188 * T**2 + 0.017998 * T**3) * a2r
    z = (2306.2181 * T + 1.09468 * T**2 + 0.018203 * T**3) * a2r
    theta = (2004.3109 * T - 0.42665 * T**2 - 0.041833 * T**3) * a2r
    ra0, dec0 = np.radians(ra_deg), np.radians(dec_deg)
    A = np.cos(dec0) * np.sin(ra0 + zeta)
    B = np.cos(theta) * np.cos(dec0) * np.cos(ra0 + zeta) - np.sin(theta) * np.sin(dec0)
    C = np.sin(theta) * np.cos(dec0) * np.cos(ra0 + zeta) + np.cos(theta) * np.sin(dec0)
    ra = (np.degrees(np.arctan2(A, B) + z)) % 360.0
    dec = np.degrees(np.arcsin(np.clip(C, -1, 1)))
    return ra, dec


def local_sidereal_deg(jd, lon_deg):
    T = (jd - 2451545.0) / 36525.0
    gmst = (280.46061837 + 360.98564736629 * (jd - 2451545.0)
            + 0.000387933 * T**2 - T**3 / 38710000.0)
    return (gmst + lon_deg) % 360.0


def alt_az(ra_deg, dec_deg, lat_deg, lst_deg):
    """ra/dec: shape (N,1); lst: shape (T,). Returns alt, az in degrees (N,T).
    Azimuth is measured from North through East. Includes refraction."""
    lat = math.radians(lat_deg)
    H = np.radians(lst_deg[None, :] - ra_deg)
    dec = np.radians(dec_deg)
    sin_alt = math.sin(lat) * np.sin(dec) + math.cos(lat) * np.cos(dec) * np.cos(H)
    alt = np.degrees(np.arcsin(np.clip(sin_alt, -1, 1)))
    az = np.degrees(np.arctan2(-np.cos(dec) * np.sin(H),
                               np.sin(dec) * math.cos(lat) - np.cos(dec) * np.cos(H) * math.sin(lat)))
    az %= 360.0
    # Bennett refraction (arcmin) - lifts objects slightly near the horizon
    h = np.maximum(alt, -1.0)
    alt = alt + np.where(alt > -1.0, 1.0 / np.tan(np.radians(h + 7.31 / (h + 4.4))) / 60.0, 0.0)
    return alt, az


def tangent_plane(ra0, dec0, ra, dec):
    """Gnomonic projection around (ra0, dec0). Returns x, y in arcmin with
    North up and East LEFT (as the sky looks), plus a mask of points in front."""
    a0, d0 = math.radians(ra0), math.radians(dec0)
    a, d = np.radians(ra), np.radians(dec)
    cosc = np.sin(d) * math.sin(d0) + np.cos(d) * math.cos(d0) * np.cos(a - a0)
    ok = cosc > 0.5
    cosc = np.where(ok, cosc, 1.0)
    xi = np.cos(d) * np.sin(a - a0) / cosc
    eta = (np.sin(d) * math.cos(d0) - np.cos(d) * math.sin(d0) * np.cos(a - a0)) / cosc
    k = 180 / math.pi * 60
    return -xi * k, -eta * k, ok


def az_in_range(az, az_from, az_to):
    """True where az lies on the clockwise arc from az_from to az_to."""
    width = (az_to - az_from) % 360.0
    if width == 0:
        width = 360.0
    return ((az - az_from) % 360.0) <= width + 1e-9


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
def load_config():
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, encoding="utf-8") as f:
            return json.load(f)
    return {}


def get_tz(name):
    if name in (None, "", "local"):
        return datetime.now().astimezone().tzinfo
    if ZoneInfo is None:
        sys.exit("Named time zones need Python 3.9+. Use e.g. --tz +02:00")
    if name[0] in "+-":
        sign = -1 if name[0] == "-" else 1
        hh, _, mm = name[1:].partition(":")
        return timezone(sign * timedelta(hours=int(hh), minutes=int(mm or 0)))
    return ZoneInfo(name)


def parse_period(args, tz):
    def parse_dt(s):
        s = s.strip().replace("T", " ")
        for fmt in ("%Y-%m-%d %H:%M", "%Y-%m-%d %H:%M:%S"):
            try:
                return datetime.strptime(s, fmt)
            except ValueError:
                pass
        return None

    day = datetime.strptime(args.date, "%Y-%m-%d").date() if args.date else date.today()

    def resolve(s):
        full = parse_dt(s)
        if full:
            return full
        t = datetime.strptime(s.strip(), "%H:%M").time()
        return datetime.combine(day, t)

    start = resolve(args.start)
    end = resolve(args.end)
    # "--from 21:00 --to 01:00" means 01:00 on the next day
    if end <= start and parse_dt(args.end) is None:
        end += timedelta(days=1)
    if end <= start:
        sys.exit("End of period must be after start.")
    return start.replace(tzinfo=tz), end.replace(tzinfo=tz)


def build_parser(cfg):
    p = argparse.ArgumentParser(
        description="Find NGC/IC/Messier objects inside a sky window during a time period.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_argument_group("observer (defaults can be stored in config.json)")
    g.add_argument("--lat", type=float, default=cfg.get("lat"), help="latitude, deg (+N)")
    g.add_argument("--lon", type=float, default=cfg.get("lon"), help="longitude, deg (+E)")
    g.add_argument("--tz", default=cfg.get("tz", "local"),
                   help="time zone, e.g. Europe/Kyiv or +02:00")
    g.add_argument("--save", action="store_true",
                   help="save location, window (--az-from/--az-to/--alt-min/--alt-max) and camera settings to config.json as new defaults")

    g = p.add_argument_group("time period (local time)")
    g.add_argument("--date", help="date of the evening, YYYY-MM-DD (default: today)")
    g.add_argument("--from", dest="start", default="21:00",
                   help="start, HH:MM or 'YYYY-MM-DD HH:MM'")
    g.add_argument("--to", dest="end", default="23:00",
                   help="end, HH:MM (next day if earlier than start) or 'YYYY-MM-DD HH:MM'")
    g.add_argument("--step", type=float, default=5, help="time step for the check, minutes")

    g = p.add_argument_group("sky window")
    g.add_argument("--az-from", default=cfg.get("az_from", "NW"),
                   help="window left edge: degrees or compass point (clockwise to --az-to)")
    g.add_argument("--az-to", default=cfg.get("az_to", "N"), help="window right edge")
    g.add_argument("--alt-min", type=float, default=cfg.get("alt_min", 10), help="lowest altitude, deg")
    g.add_argument("--alt-max", type=float, default=cfg.get("alt_max", 45), help="highest altitude, deg")
    g.add_argument("--whole", action="store_true",
                   help="only objects that stay in the window for the WHOLE period")
    g.add_argument("--min-minutes", type=float, default=0,
                   help="minimum time the object spends in the window")

    g = p.add_argument_group("object filters")
    g.add_argument("--catalogs", default="ngc,ic,m",
                   help="any of: ngc, ic, m (m = all Messier objects)")
    g.add_argument("--messier-only", action="store_true", help="shortcut for --catalogs m")
    g.add_argument("--types", default="galaxy,cluster,nebula",
                   help="any of: galaxy, cluster, nebula, other")
    g.add_argument("--max-mag", type=float, default=12.0, help="faintest magnitude")
    g.add_argument("--max-sb", type=float, help="faintest surface brightness, mag/arcsec^2")
    g.add_argument("--min-size", type=float, default=0, help="smallest major axis, arcmin")
    g.add_argument("--no-mag", action="store_true",
                   help="also include objects that have no magnitude in the catalog")
    g.add_argument("--favorites", default=cfg.get("favorites", FAVORITES_FILE),
                   help="file with favorite objects, one per line (e.g. 'M 101', 'NGC 6946')")
    g.add_argument("--favorites-only", action="store_true",
                   help="only list favorites (other filters like --max-mag are ignored for them)")

    g = p.add_argument_group("camera / telescope (for frame size and the HTML frame view)")
    g.add_argument("--focal", type=float, default=cfg.get("focal", 750),
                   help="telescope focal length, mm (150 mm f/5 = 750)")
    g.add_argument("--pixel", type=float, default=cfg.get("pixel", 2.9), help="pixel size, microns")
    g.add_argument("--sensor", default=cfg.get("sensor", "1936x1096"), help="sensor size, pixels WxH")
    g.add_argument("--camera-name", default=cfg.get("camera_name", "ASI290MC + 150/750"),
                   help="label shown on the frame")

    g = p.add_argument_group("output")
    g.add_argument("--html", help="write an interactive HTML report with the camera-frame view")
    g.add_argument("--sort", default="mag",
                   choices=["mag", "sb", "size", "name", "alt", "time", "enter"],
                   help="sort order (alt = altitude at start, time = minutes in window)")
    g.add_argument("--limit", type=int, default=0, help="show only first N rows (0 = all)")
    g.add_argument("--csv", help="also write results to this CSV file")
    return p


def fmt(v, spec, empty="-"):
    return empty if v is None else format(v, spec)


def main():
    cfg = load_config()
    args = build_parser(cfg).parse_args()

    if args.lat is None or args.lon is None:
        sys.exit("Observer location needed: --lat 50.45 --lon 30.52 "
                 "(add --save to remember it in config.json)")
    if args.save:
        cfg.update({"lat": args.lat, "lon": args.lon, "tz": args.tz,
                    "az_from": args.az_from, "az_to": args.az_to,
                    "alt_min": args.alt_min, "alt_max": args.alt_max,
                    "focal": args.focal, "pixel": args.pixel, "sensor": args.sensor,
                    "camera_name": args.camera_name})
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)
        print(f"Saved settings to {CONFIG_FILE}")

    tz = get_tz(args.tz)
    start, end = parse_period(args, tz)
    az_from, az_to = parse_az(args.az_from), parse_az(args.az_to)

    # --- select candidate objects -----------------------------------------
    cats = {c.strip().lower() for c in args.catalogs.split(",")}
    if args.messier_only:
        cats = {"m"}
    types = {t.strip().lower().rstrip("s") for t in args.types.split(",")}

    sw, sh = (int(v) for v in args.sensor.lower().split("x"))
    scale = 206.265 * args.pixel / args.focal              # arcsec per pixel
    frame_w, frame_h = sw * scale / 60, sh * scale / 60     # arcmin

    catalog = load_catalog()
    favs = load_favorites(args.favorites)
    if args.favorites_only and not favs:
        sys.exit(f"No favorites found in {args.favorites}. Add one object per line, e.g. 'M 101'.")

    objs = []
    for o in catalog:
        o["fav"] = is_favorite(o, favs)
        if args.favorites_only:
            if o["fav"]:
                objs.append(o)
            continue
        in_cat = (o["cat"] in cats) or ("m" in cats and o["messier"])
        if not in_cat or o["group"] not in types:
            continue
        if o["mag"] is None:
            if not args.no_mag:
                continue
        elif o["mag"] > args.max_mag:
            continue
        if args.max_sb is not None and (o["sb"] is None or o["sb"] > args.max_sb):
            continue
        if args.min_size and (o["maj"] or 0) < args.min_size:
            continue
        objs.append(o)
    if not objs:
        sys.exit("No catalog objects match the filters.")

    # --- compute positions over the period --------------------------------
    n_steps = max(1, int(math.ceil((end - start).total_seconds() / 60 / args.step)))
    times = [start + (end - start) * i / n_steps for i in range(n_steps + 1)]
    jds = np.array([julian_date(t.astimezone(timezone.utc)) for t in times])
    lst = local_sidereal_deg(jds, args.lon)

    ra = np.array([o["ra"] for o in objs])
    dec = np.array([o["dec"] for o in objs])
    ra, dec = precess(ra, dec, jds[len(jds) // 2])
    alt, az = alt_az(ra[:, None], dec[:, None], args.lat, lst)

    inside = (alt >= args.alt_min) & (alt <= args.alt_max) & az_in_range(az, az_from, az_to)
    step_min = (end - start).total_seconds() / 60 / n_steps

    # sample indices for the sky path drawn in the HTML report (max ~60 points)
    k = max(1, math.ceil(len(times) / 60))
    track_idx = list(range(0, len(times), k))
    if track_idx[-1] != len(times) - 1:
        track_idx.append(len(times) - 1)

    results = []
    for i, o in enumerate(objs):
        row = inside[i]
        if not row.any():
            continue
        if args.whole and not row.all():
            continue
        idx = np.nonzero(row)[0]
        # minutes in window: count samples, each worth one step (endpoints half)
        weights = np.full(len(row), step_min)
        weights[0] = weights[-1] = step_min / 2
        minutes = float(weights[row].sum())
        if minutes < args.min_minutes:
            continue
        results.append({
            **o,
            "alt_start": float(alt[i, 0]), "alt_end": float(alt[i, -1]),
            "az_start": float(az[i, 0]), "az_end": float(az[i, -1]),
            "enter": times[idx[0]], "leave": times[idx[-1]],
            "minutes": minutes,
            "track": [[round(float(alt[i, j]), 1), round(float(az[i, j]), 1)] for j in track_idx],
        })

    # --- other catalog objects near each result (for the frame view) ------
    all_ra = np.array([c["ra"] for c in catalog])
    all_dec = np.array([c["dec"] for c in catalog])
    for r in results:
        # search far enough for the widest zoom of the HTML view
        view = max(frame_w, (r["maj"] or 0) * 0.8) * 1.1 * 6
        radius = view * 0.62 + 30
        x, y, ok = tangent_plane(r["ra"], r["dec"], all_ra, all_dec)
        near = np.nonzero(ok & (np.hypot(x, y) <= radius))[0]
        nb = []
        for j in near:
            c = catalog[j]
            if c["name"] == r["name"] or (c["messier"] and c["messier"] == r["messier"]):
                continue
            nb.append({"name": c["name"] + (f" (M {c['messier']})" if c["messier"] and c["cat"] != "other" else ""),
                       "type": c["type"], "group": c["group"], "code": c["code"],
                       "mag": c["mag"], "maj": c["maj"], "min": c["min"], "pa": c["pa"],
                       "x": round(float(x[j]), 3), "y": round(float(y[j]), 3)})
        r["nb"] = nb
        r["n_in_frame"] = sum(1 for n in nb if abs(n["x"]) <= frame_w / 2 and abs(n["y"]) <= frame_h / 2)

    big = 99.0
    sort_keys = {
        "mag": lambda r: (r["mag"] if r["mag"] is not None else big),
        "sb": lambda r: (r["sb"] if r["sb"] is not None else big),
        "size": lambda r: -(r["maj"] or 0),
        "name": lambda r: r["name"],
        "alt": lambda r: -r["alt_start"],
        "time": lambda r: -r["minutes"],
        "enter": lambda r: r["enter"],
    }
    results.sort(key=sort_keys[args.sort])
    if args.limit:
        shown = results[:args.limit]
    else:
        shown = results

    # --- print -------------------------------------------------------------
    def compass(a):
        pts = list(COMPASS)
        return pts[int((a + 11.25) // 22.5) % 16]

    print()
    print(f"Observer: lat {args.lat:+.3f}, lon {args.lon:+.3f}   Time zone: {start.tzname()}")
    print(f"Period:   {start:%Y-%m-%d %H:%M} -> {end:%Y-%m-%d %H:%M}  ({(end - start).total_seconds() / 3600:.1f} h)")
    print(f"Window:   azimuth {az_from:.0f}..{az_to:.0f} deg ({compass(az_from)}..{compass(az_to)}), "
          f"altitude {args.alt_min:.0f}..{args.alt_max:.0f} deg"
          + ("   [whole period]" if args.whole else ""))
    if args.favorites_only:
        print(f"Filters:  favorites only ({len(favs)} in {os.path.basename(args.favorites)})")
    else:
        print(f"Filters:  catalogs={','.join(sorted(cats))}  types={','.join(sorted(types))}  "
              f"mag<={args.max_mag}" + (f"  SB<={args.max_sb}" if args.max_sb else "")
              + (f"  size>={args.min_size}'" if args.min_size else ""))
    print(f"Camera:   {args.camera_name}: frame {frame_w:.1f}' x {frame_h:.1f}', {scale:.2f}\"/px")
    print(f"Found:    {len(results)} objects" + (f" (showing {len(shown)})" if len(shown) < len(results) else ""))
    print()

    headers = ["Object", "Also", "Type", "Mag", "SB", "Size '", "Fill", "+Frm", "Alt start", "Alt end",
               "Az start", "Az end", "In window", "Min", "Con", "Name"]
    rows = []
    for r in shown:
        size = "-" if r["maj"] is None else (
            f"{r['maj']:.1f}" if r["min"] == r["maj"] else f"{r['maj']:.1f}x{r['min']:.1f}")
        rows.append([
            ("★ " if r["fav"] else "") + r["name"],
            ", ".join(r["aliases"]),
            r["type"],
            fmt(r["mag"], ".1f") + (r["mag_band"] if r["mag_band"] == "B" else ""),
            ("~" if r["sb_est"] else "") + fmt(r["sb"], ".1f"),
            size,
            "-" if r["maj"] is None else f"{100 * r['maj'] / frame_w:.0f}%",
            f"+{r['n_in_frame']}" if r["n_in_frame"] else "",
            f"{r['alt_start']:5.1f}",
            f"{r['alt_end']:5.1f}",
            f"{r['az_start']:5.1f}",
            f"{r['az_end']:5.1f}",
            f"{r['enter']:%H:%M}-{r['leave']:%H:%M}",
            f"{r['minutes']:.0f}",
            r["const"],
            r["common"][:28],
        ])
    if rows:
        widths = [max(len(h), *(len(x[c]) for x in rows)) for c, h in enumerate(headers)]
        right = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13}
        line = lambda cells: "  ".join(
            (c.rjust(w) if i in right else c.ljust(w)) for i, (c, w) in enumerate(zip(cells, widths)))
        print(line(headers).rstrip())
        print("  ".join("-" * w for w in widths))
        for x in rows:
            print(line(x).rstrip())
        print()
        print("Mag: V magnitude (B = blue magnitude when V is unknown).  "
              "SB: surface brightness, mag/arcsec^2 (~ = estimated from mag and size).")
        print("Alt/Az in degrees at the start and end of the period (azimuth from North, clockwise).  "
              "In window: first-last time inside the window.")
        print("Fill: major axis as % of frame width.  +Frm: other NGC/IC objects inside the frame "
              "(camera long side East-West).")

    if args.html:
        from report_html import write_report
        rows = []
        for i, r in enumerate(results):
            rows.append({
                "i": i, "name": r["name"], "aliases": r["aliases"], "type": r["type"],
                "group": r["group"], "code": r["code"], "mag": r["mag"], "mag_band": r["mag_band"],
                "sb": None if r["sb"] is None else round(r["sb"], 2), "sb_est": r["sb_est"],
                "maj": r["maj"], "min": r["min"], "pa": r["pa"],
                "alt_start": round(r["alt_start"], 1), "alt_end": round(r["alt_end"], 1),
                "az_start": round(r["az_start"], 1), "az_end": round(r["az_end"], 1),
                "enter": f"{r['enter']:%H:%M}", "leave": f"{r['leave']:%H:%M}",
                "enter_ts": r["enter"].timestamp(), "minutes": r["minutes"],
                "const": r["const"], "common": r["common"],
                "ra": round(r["ra"], 5), "dec": round(r["dec"], 5),
                "ra_s": r["ra_s"], "dec_s": r["dec_s"],
                "nb": r["nb"], "n_in_frame": r["n_in_frame"],
                "track": r["track"],
            })
        write_report(args.html, {
            "period": f"{start:%Y-%m-%d %H:%M} → {end:%H:%M} ({start.tzname()})",
            "window": f"az {az_from:.0f}–{az_to:.0f}° ({compass(az_from)}–{compass(az_to)}), "
                      f"alt {args.alt_min:.0f}–{args.alt_max:.0f}°",
            "lat": args.lat, "lon": args.lon,
            "win": {"az_from": az_from, "az_to": az_to, "alt_min": args.alt_min, "alt_max": args.alt_max},
            "track_times": [f"{times[j]:%H:%M}" for j in track_idx],
            "camera": {"name": args.camera_name, "fw": frame_w, "fh": frame_h, "scale": scale},
            "rows": rows,
            "favorites": [r["name"] for r in results if r["fav"]],
            "constellations": constellation_payload(
                jds[0], lst[0], args.lat, (az_from + ((az_to - az_from) % 360.0 or 360.0) / 2) % 360),
        })
        print(f"\nWrote HTML report to {args.html}")

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["object", "also", "group", "type", "mag", "mag_band", "surface_brightness",
                        "sb_estimated", "size_major_arcmin", "size_minor_arcmin",
                        "alt_start", "alt_end", "az_start", "az_end",
                        "window_enter", "window_leave", "minutes_in_window",
                        "constellation", "common_name", "ra_j2000_deg", "dec_j2000_deg"])
            for r in results:
                w.writerow([r["name"], "; ".join(r["aliases"]), r["group"], r["type"],
                            fmt(r["mag"], ".2f", ""), r["mag_band"], fmt(r["sb"], ".2f", ""),
                            "yes" if r["sb_est"] else "", fmt(r["maj"], ".2f", ""),
                            fmt(r["min"], ".2f", ""), f"{r['alt_start']:.1f}",
                            f"{r['alt_end']:.1f}", f"{r['az_start']:.1f}", f"{r['az_end']:.1f}",
                            f"{r['enter']:%Y-%m-%d %H:%M}", f"{r['leave']:%Y-%m-%d %H:%M}",
                            f"{r['minutes']:.0f}", r["const"], r["common"],
                            f"{r['ra']:.5f}", f"{r['dec']:.5f}"])
        print(f"\nWrote {len(results)} rows to {args.csv}")


if __name__ == "__main__":
    main()
