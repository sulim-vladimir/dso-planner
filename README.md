# DSO Window Planner

Finds NGC / IC / Messier objects that pass through a patch of sky (by default
**azimuth NW–N, altitude 10–45°**) during a time period, and shows how each one
fits in the camera frame (default: **ASI290MC on a 150 mm f/5 reflector**).

![HTML report: object list and M51 in the ASI290MC frame](docs/screenshot.png)

Needs Python 3.9+ and numpy. Catalog: [OpenNGC](https://github.com/mattiaverga/OpenNGC) (CC-BY-SA-4.0), in `data/`.
Constellation lines and names: [Stellarium](https://stellarium.org) (western sky culture);
star positions: [HYG database](https://github.com/astronexus/HYG-Database) (CC BY-SA 4.0).

## Quick start

```sh
# first time: give your location once and save it (+N latitude, +E longitude)
./dso_window.py --lat 50.45 --lon 30.52 --tz Europe/Kyiv --save

# change the window (e.g. maximum altitude 35°) and keep it as the new default
./dso_window.py --alt-max 35 --save

# tonight 21:00 → 01:00, terminal table + interactive HTML report
./dso_window.py --from 21:00 --to 01:00 --html tonight.html

# a specific date, galaxies only, brighter than mag 11, longest in the window first
./dso_window.py --date 2026-11-14 --from 20:30 --to 23:00 --types galaxy --max-mag 11 --sort time
```

Open the HTML file in a browser and click any row to see it in the camera frame.

## What you get

| Column | Meaning |
|---|---|
| Type | galaxy / cluster / nebula (with subtype, e.g. *Galaxy SABc*, *Planetary nebula*) |
| Mag | V magnitude (`B` suffix = blue magnitude when V is unknown) |
| SB | surface brightness, mag/arcsec² (`~` = estimated from mag and size) |
| Size | major × minor axis, arcmin |
| Fill | major axis as % of the frame width |
| +Frm | other NGC/IC objects that also land in the frame |
| Alt/Az start, end | position at the start and end of the period |
| In window / Min | first–last time inside the window, and minutes spent there |

## HTML frame view

- Schematic of the object (true size and orientation) inside the 25.7′ × 14.6′ frame, with
  every other NGC/IC object nearby — galaxies, clusters, nebulae.
- **Rotation** slider = camera rotation; **Zoom** slider = field around the frame.
- **Where in the sky**: a view facing the window (stereographic projection, as in Stellarium, so shapes are true) with the horizon, altitude lines, the window
  outlined, the object's position at the start of the period and its path (hour marks) until the end.
  Faint dots are the other objects of the list at the start; click one to select it.
  Constellation lines (at the start of the period) help to find the spot; hover a star for its name.
- Optional photo underlay to compare with the real sky; links to Aladin, SIMBAD.
  *Auto* uses Legacy Survey (color, sharp) and switches to DSS2 red (NASA SkyView, whole sky)
  for Milky Way objects, which Legacy Survey does not cover. The source can also be picked by hand.
- **Night mode** turns everything red.
- **↑ / ↓** (or `j` / `k`) move to the previous / next object in the table.
- **Copy** button next to the object name (or press `c`) copies its catalog name, e.g. `NGC 6229`,
  for pasting into a mount / planetarium app.
- **☆ / ★ favorites**: click the star in a row or next to the name (or press `f`);
  the **★ Favorites** chip shows only favorites. See [Favorites](#favorites).
- Open a given object directly: `tonight.html#M 101`.

<img src="docs/frame-photo.png" alt="NGC 5866 in the frame with the survey photo underlay" width="420">

*NGC 5866 with the photo underlay on: the schematic outline matches the real galaxy.*

The page layout lives in `templates/report.html`.

## Favorites

- **In the report**, stars are remembered by your browser, so they carry over to new reports
  (tested in Chrome; they are kept per browser, not shared between computers).
- **In the terminal**, list favorites in `favorites.txt` next to the script, one per line —
  any of `M 101`, `NGC 6946`, `ngc6946`, or a common name like `Cat's Eye Nebula`; `#` starts a comment.
  They get a ★ in the table and are starred in the report too.
- **Copy list** in the report copies your starred objects, one per line, ready to paste into `favorites.txt`.

```sh
./dso_window.py --from 21:00 --to 01:00 --favorites-only       # only favorites in the window tonight
```

## Useful options

| Option | Default | |
|---|---|---|
| `--az-from`, `--az-to` | `NW`, `N` | window edges, degrees or compass points, clockwise |
| `--alt-min`, `--alt-max` | 10, 45 | window altitude limits |
| `--whole` | off | only objects that stay in the window the whole period |
| `--min-minutes` | 0 | minimum time spent in the window |
| `--catalogs` / `--messier-only` | `ngc,ic,m` | catalogs to search |
| `--types` | `galaxy,cluster,nebula` | add `other` for stars etc. |
| `--max-mag`, `--max-sb`, `--min-size` | 12, –, 0 | brightness / size filters |
| `--no-mag` | off | also list objects with no catalog magnitude (many nebulae) |
| `--favorites-only` | off | only objects listed in `favorites.txt` (ignores mag/type filters) |
| `--favorites FILE` | `favorites.txt` | use a different favorites file |
| `--sort` | `mag` | `mag sb size name alt time enter` |
| `--csv FILE` | | save results as CSV |
| `--focal`, `--pixel`, `--sensor` | 750, 2.9, 1936x1096 | telescope/camera (can be saved in `config.json`) |

Run `./dso_window.py -h` for the full list. `./update_catalog.sh` refreshes the catalog and the constellation star positions.
