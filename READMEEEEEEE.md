# Transit Hunter — Physics Club Team

**HACK4DEV IRAQ 2026 — Exoplanet Data Challenge**
**Challenge A — Transit Hunter**

---

## The scientific question

**Can a 6-inch robotic telescope recover the known transit of an exoplanet, and what is the precision limit of our pipeline?**

---

## Headline result

We did not obtain a confirmed transit detection, and we state that up front.

Across the campaign, the photometric scatter we could achieve (1.08% in our best night, 4–7% typically) is comparable to or larger than the published transit depths of these planets (1.5–2.9%). We therefore report **non-detections with quantified sensitivity**, plus one marginal 3.2σ dip on Qatar-1 (2026-08-21) that our own reverse test refuses to certify.

The pipeline, the diagnostics, and the rejection logic are the deliverable.

---

## Why we chose this challenge

We chose Transit Hunter because it engaged us most directly with the physics and astronomy we care about. The other challenges could be approached as data problems; this one required thinking about starlight, atmosphere, and orbital motion, and turning raw telescope images into a physical measurement.

The data pointed the same way. We were given raw FITS frames with no labels and no feature table, so the machine-learning challenges could not be applied directly. What the data does contain is exactly what a transit search needs: a dense time series of the same field across a night.

And because the planets here are already confirmed and published, we had an independent reference to check ourselves against — which mattered more to us than the chance of finding something new.

---

## The problem we solved

If we look at the target star alone and see its brightness drop, we cannot be sure. It might be a planet passing in front of it, or it might be a flaw in how we measured. The idea we used is simple: a planet passes in front of one star, and it cannot make any star brighter. It blocks light; it does not add light. So we applied the same analysis to five other stars in the image and looked at what happened to them at the same moment. If the target alone dipped while the rest stayed flat, the signal is real. If other stars dipped by the same amount, and some stars grew brighter, then the problem is in the measurement, not in the sky.

This actually saved us. On one night we found a dip of 2.08%, and the published value for that planet is 2.20%, a near-perfect match, and we were about to report it as a detection. Then we ran the test and found another star in the same image that had grown brighter by more than our target had dimmed. Since a planet cannot do that, we knew our method produces numbers like these on its own, and that the match was a coincidence. So we did not report a detection. This test is what kept us from presenting a false result.

---

## The most important problem in the data

The star field moves during the night. The telescope does not hold perfectly still, and we measured drifts of up to 75 pixels. This matters because we measure brightness by drawing a small circle around the star and adding up the light inside it. If the circle stays fixed while the star drifts away, by mid-night we are measuring empty sky — and the brightness appears to collapse, which looks exactly like a transit.

The obvious fix is to search for the brightest pixel near the star's last position. We tried it and it failed: our target is faint, so when a brighter neighbour drifted into the search box the tracking jumped onto it, and the measured brightness leapt from 782 to 44,602 in one frame.

What worked was following the whole field instead of one star. We detect all the stars in each frame and match their geometric pattern against a reference — a triangle of three stars has a shape that does not repeat by chance, so a single confused star gets corrected by the group. Once we know how far the field shifted, we move every measurement circle by the same amount. This kept all 74 frames of that night usable.

---

## Dataset

| Item | Value |
|---|---|
| Source | MicroObservatory Image Directory (Harvard/Smithsonian) |
| Telescope | `Cecilia`, 6-inch robotic |
| Total FITS files | 1741 (1681 science + 60 darks) |
| Targets | 8 |
| Observing sessions | 22 |
| Date range | 2026-08-06 → 2026-09-05 |
| Image properties | 650 × 500 px, 5.0 arcsec/px, 60 s exposure, Clear filter, 12-bit (max 4095) |

**No flat fields exist in the package.** Pixel-to-pixel sensitivity is uncorrected — a stated limitation of every result below.

**No WCS in the headers.** The `RA`/`DEC` keywords give the telescope pointing, not the target position, and the image orientation is undocumented. Solving this shaped the entire project.

---

## Tools

| Tool | Use |
|---|---|
| Python 3 | All analysis |
| `astropy` | FITS I/O, headers, sigma-clipped statistics |
| `photutils` | Star detection (`DAOStarFinder`), aperture photometry |
| `numpy` / `pandas` | Numerical work, tables, logs |
| `matplotlib` | All figures |
| Gaia DR3 catalog | Plate-solving the fields and selecting comparison stars |
| NASA Exoplanet Archive / ETD | Published transit depths and ephemerides (external source, disclosed) |

---

## How to run

### Requirements

```bash
pip install numpy scipy pandas matplotlib astropy photutils
```

### Expected data layout

```
OneStar/
├── database/
│   ├── observations/YYYY-MM-DD/TARGET/session_01/*.fits   # 1681 science frames
│   └── calibration/YYYY-MM-DD/*.fits                      # 60 dark frames
└── results/                                               # all outputs land here
```

The raw archive (~1.1 GB) is not in the repository. Re-download from the MicroObservatory Image Directory (ExoPlanets, last 30 days), or request it from the team.

### Run order

```bash
python database/index.py       # build the dataset index (1741 rows)
python find_real_targets.py    # plate-solve every frame against Gaia   [long]
python check_cutouts_real.py   # verify a real star sits at the solved position
python get_star_positions.py   # sub-pixel centroid refinement
python check_transits.py       # which sessions actually contain a transit
python photometry.py           # differential photometry, per-session light curves
python combine_star.py         # per-target campaign reports
```

Standalone single-night pipeline (independent of the Gaia solve):

```bash
python hassan                  # edit BASE / NIGHT / TARGET at the top first
```

Frame cleaning and clean-subset construction:

```bash
python clean_frames.py         # score all 1681 frames on image content
python build_clean_dataset.py  # copy the 5 selected nights into database_clean/
```

---

## Scientific interpretation

**The field can be solved while the target star is still undetectable.** Gaia pattern-matching solved 644 of 1681 frames and recovered the undocumented image orientation (`flip-y`, 585 votes against 51 and 8). But stacking all frames of a session at the solved target position showed a real star in only 7 of 17 solved sessions. Per frame, a detected star coincides with the solved target position in just 1 of 1681 frames. All photometry below is therefore performed at the *solved catalog* position, not at an independently detected star — and we say so rather than hiding it.

**Every non-detection is a non-detection of an event that should have been there.** Propagating published ephemerides to each observing window showed that all 22 sessions overlap a predicted transit, and 21 of 22 contain it fully. MicroObservatory schedules these targets on transit nights. That is what makes our sensitivity numbers meaningful: we were not looking in the wrong place.

**The noise floor is real, the target noise is worse.** Check-star scatter — a comparison star analysed as if it were the target, our control experiment — reaches 0.3–3.1%, and in the best sessions falls *below* the expected transit depth. But target scatter is consistently several times worse. That gap is the signature of the host star sitting at or below the detection limit at 5 arcsec per pixel.

**Our best night, quantified.** Qatar-1 on 2026-08-21: 74 frames, all passing quality cuts, 3.62 h span, airmass 1.20 → 1.33, field drift +39 px in x and +67 px in y. Final scatter after detrending 1.08%. Measured depth +0.74% at 3.22σ. Published depth for Qatar-1b is 2.0%, giving a sensitivity ratio of 1.9:1 — below the 2:1 our pipeline requires before it will certify anything.

**The reverse test decided it.** Applying the identical analysis to six other stars in the same frames, the target dimmed by 0.648% while a control star dimmed by 0.368% — only 1.8× separation, and two other stars *brightened*. A planet cannot brighten a star. We report the dip as a marginal, uncertified candidate, not a transit.

**What the result means.** With 60-second exposures from a 6-inch telescope at 5 arcsec per pixel, on host stars near the detection limit, and with no flat-field correction available, a 2% signal cannot be separated from residual systematics. This is a measured ceiling, not an assumed one, and reporting it honestly is the scientific content of our submission.

---

## Limitations

1. **No flat fields** in the MicroObservatory distribution. Partially mitigated by detrending relative flux against (x, y) position — an approximation, not a calibration.
2. **The target is measured at the solved catalog position**, not at an independently detected star, in almost all frames.
3. **Ephemeris propagation error:** literature T0 values propagated over a decade give ~10–30 min uncertainty in predicted mid-transit; BJD_TDB versus UTC adds up to ~8 min. Adequate for yes/no window overlap, not for timing measurements — we make none.
4. **5 arcsec per pixel is coarse.** Neighbouring stars blend into the apertures; contamination dilutes any real depth.
5. **Orientation was determined statistically**, not from documentation. We did not force a single convention on the data.
6. **62% of frames are unsolved**, mostly weather. One entire night (TRES-1, 2026-08-06) has `WEATHER = 0` and no detectable stars at all.

---

## What we would do next

- **Bin in time** to 10–15 minute cadence before searching for the dip. Binned check-star scatter already reaches sub-1% in several sessions and would push sensitivity below the 2% depths.
- **Weight comparison stars** by inverse variance rather than summing equally, and require at least 4 comparisons instead of 2.
- **Fit a proper transit model** (Mandel–Agol with limb darkening) to the binned curve with the ephemeris fixed, and report a depth *upper limit* per night.
- **Build a master flat** from the science frames themselves, using the median of many dithered frames of different fields.
- **Inject synthetic transits** of known depth into the data and measure how often the pipeline recovers them — turning our sensitivity estimate into a measurement.

---

## Reading the campaign figures

Three campaign reports are reproduced here — Qatar-1, TRES-3 and TRES-5. Each has the same structure, and each is read the same way. Because this is where the actual evidence lives, we explain them in full.

### How to read the layout

Every figure has two levels.

**The top strip is the campaign overview.** The horizontal axis is real calendar date across the whole observing run, 6 August to 5 September. Each marker is one night, placed at its true date, and its vertical position is the median relative flux of that night. A filled black square means the night produced usable photometry; a hollow square means it was measured but flagged noisy; a grey cross labelled *no data* marks a night that was scheduled and observed but produced no measurable light curve. The vertical error bar is the spread of that night's measurements.

**The lower panels are the individual nights.** The horizontal axis is UT clock time, not frame number, so gaps in the data appear as real gaps. Black points with error bars are the binned light curve; faint grey points behind them are the individual frames. The orange line is the **predicted transit model** — a trapezoid built from the published depth and duration and placed at the ephemeris-predicted time. It is drawn on every panel whether or not data exist there, so the reader can see immediately whether our observations even covered the event.

**One critical caveat in the titles.** Each figure states how its nights were normalized. TRES-5 reads *4 shared comps* and TRES-3 reads *2 shared comps*, meaning those nights share comparison stars and sit on a common brightness scale — so comparing them to each other is valid. Qatar-1 reads *nights normalized individually (no shared comps)*: no comparison star was usable across all three nights, so each night was normalized to its own median. **This makes the top strip of the Qatar-1 figure non-comparable between nights.** The value 1.00 there is imposed by normalization, not measured. We label it rather than hiding it.

### TRES-5 — our best data

Two nights produced light curves: 26 August with 6 frames, and 29 August with 62.

The overview strip shows 08-26 sitting at 0.957 and 08-29 at 1.004 — a 4.7% difference between two nights of the same star. A star does not change brightness by 5% between nights, so this difference is not astrophysical. It is an artefact of 08-26 surviving with only six frames, all of them clustered at the tail end of the observation after the predicted transit had already ended. With that few points, the night median is meaningless as a brightness measurement.

The 08-26 panel makes the timing failure visible: the orange transit model runs from roughly 07:40 to 09:10 UT, and our six surviving points begin at 09:30. **The event was over before our usable data starts.** No conclusion about a transit is possible from this night, and we draw none.

The 08-29 panel is the only one in the entire campaign where the data genuinely covers the predicted window. Sixty-two frames span 06:00 to 09:30, and the orange model sits at roughly 07:00 to 08:30, fully inside the observation. There is a visible group of low points around 07:45 to 08:15 reaching 0.94, which is close to the predicted depth. But reading the rest of the panel settles it: points before the window reach 1.09, and a point shortly after the window drops to 0.87. The overall scatter of the night is far larger than the dip, so the low group does not stand out from the noise — it sits inside it. Per-point error bars here are 3 to 5%, against a predicted depth of 2.2%. **The uncertainty on a single measurement exceeds the signal being sought.**

### Qatar-1 — measured, but not comparable

Three nights were scheduled. One produced nothing at all.

**18 August** is flagged NOISY by the pipeline, and the panel shows why: error bars extend past 20%, several points run off the top and bottom of the plot, and grey individual frames are scattered from 0.87 to 1.13. The predicted transit sits around 10:00–10:30 UT and there are points near it, but at this noise level the panel carries no information about a 2% dip.

**21 August** is the cleanest night of the whole campaign. 39 binned points from 04:40 to 08:20 UT, scatter around ±4%, and a visible group of lower points between 05:40 and 07:00 reaching 0.97. The predicted transit window sits inside that region. This is the 3.22σ, 0.74% candidate quantified in the results table above — and it is also the one our reverse test refused to certify, because a control star in the same frames moved by more than half as much.

**28 August** produced no measurable light curve. The panel carries the explicit reason: *clouds / unsolved field*. This is consistent with our independent frame-cleaning run, which measured a median sky background of 2250 counts on this night against a typical 410 — roughly five times normal, the signature of heavy cloud or bright moonlight. The night was correctly identified as unusable by two separate methods.

### TRES-3 — three empty nights

Four nights were scheduled. Three produced no data at all — 10 August, 14 August and 31 August all carry the *clouds / unsolved field* label, meaning the star field could not be matched against the catalog in enough frames to build a curve.

The single surviving night, 27 August, produced **14 binned points**. The panel shows the decisive problem clearly: the orange transit model runs from about 06:45 to 08:00 UT, and every one of our fourteen points falls after 08:15. **The transit had ended before our usable data began.** The scatter in those points, roughly ±5%, is in any case larger than the 2.5% depth we would be looking for.

We include this figure precisely because it produced nothing. It documents a scheduled observation that could not answer the question, and the reason is visible rather than asserted.

### What the figures establish as a set

Read together, the three reports show three independent reasons why this campaign could not confirm a transit, and each is visible on the page rather than claimed in text.

**Coverage.** Of the nights that produced any light curve at all, only one — TRES-5 on 29 August — actually overlaps the predicted transit window with usable data. The orange model drawn on every panel is what makes this checkable at a glance.

**Attrition.** Nights that passed image-quality screening still collapsed at the photometric stage. TRES-5 on 26 August entered with 74 accepted frames and produced six; three of four TRES-3 nights produced none. Image quality and geometric registration are distinct failure modes, and the figures show the second one operating after the first had passed.

**Precision.** Per-point error bars run 3 to 5% across all three targets, against predicted depths of 2.0 to 2.5%. The error on a single measurement is larger than the thing being measured. No arrangement of these points can establish a detection, and the figures show that limit directly.

This is why we report a measured sensitivity ceiling rather than a detection. The figures are the evidence for that claim, not an illustration of it.

---

## Use of AI tools

We used an AI assistant during this project, and we disclose exactly how.

**Where we used it.** As a coding assistant while writing and debugging the Python pipeline; as a tutor for concepts none of us had worked with before (aperture photometry, sky-background subtraction via an annulus, airmass and atmospheric extinction, differential photometry, statistical significance); and as a language editor for this README and the written explanations.

**Where we did not use it.** No AI model is part of the solution itself. The pipeline contains no machine learning of any kind — it is classical aperture photometry using `astropy` and `photutils`, both open-source and cited above. No result, number, or figure in this submission was produced by an AI model. Every value reported here came from running our own code on the data.

**How we verified it.** Every piece of code we kept, we ran ourselves and checked against the data. Several suggestions turned out to be wrong and we discarded them after testing:

- An early approach identified the target as the brightest star in the field. We tested it and found that star sat 235 pixels from the field centre and was saturated in two-thirds of the frames. We replaced the rule with "closest unsaturated star to the field centre."
- A local brightest-pixel tracker was suggested for following the star across the night. We ran it and it jumped onto a neighbouring star, producing a flux leap from 782 to 44,602. We replaced it with the three-star geometric pattern matching described above.
- An early assessment concluded there was no transit anywhere in the package. We checked it against published ephemerides and found the opposite: all 22 sessions were scheduled on transit nights. The assessment was wrong and we corrected it.

**What we understand.** Every stage of the pipeline can be explained by a member of the team: why the sky annulus uses a median rather than a mean, why the airmass detrend must be fitted outside the transit window, why the reverse test is decisive, and why we report a non-detection rather than a 3.2σ candidate. We did not keep any code we could not explain.

**External resources disclosed.** Published transit depths and ephemerides from the NASA Exoplanet Archive and the Exoplanet Transit Database; field star catalogs from Gaia DR3; open-source libraries `astropy`, `photutils`, `numpy`, `pandas`, `scipy`, `matplotlib`.

---

## A note on what we did not claim

We present no feature in our results as a detected transit. The planets in this package are confirmed and published; our aim was to test whether our pipeline could recover a known signal and to measure the limits of its sensitivity. Where our own tests refused to certify a candidate, we report the refusal.

External data sources used: NASA Exoplanet Archive / Exoplanet Transit Database (published depths and ephemerides), Gaia DR3 (field star catalogs).

---

**Physics Club Team** — HACK4DEV IRAQ 2026
