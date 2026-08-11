# Study Guide — Understand This Project Well Enough to Defend It

**Document ID:** SG-001
**Audience:** you — someone who can run and interpret simulations but doesn't want
to derive the math from scratch. This guide gets you to the point where you can
explain *what the project does, why, and what it found*, in your own words, and
handle the follow-up questions honestly.

Read this top to bottom once, then skim the "Say it like this" boxes before the
interview.

---

## 1. The one-paragraph version (memorize the shape of this)

> "I modeled the temperature inside the fuel block of a helium-cooled nuclear
> microreactor — the kind Radiant builds. The block makes heat in its fuel and has
> to get rid of it. I computed the temperature everywhere for two cases: normal
> operation with the coolant flowing, and a loss-of-cooling accident where the only
> way heat escapes is by radiating off the outer surface. The key result is that
> even with all forced cooling gone, the fuel stays about 1240 degrees below its
> safety limit — it's passively safe. I checked the model against problems with
> known answers to make sure the numbers are trustworthy."

Everything below is just backup for that paragraph.

---

## 2. The physics, in plain language

Four things happen to heat in this block:

1. **Generation** — the fuel compacts produce heat from fission (a "volumetric heat
   source," like the whole material glowing warm from the inside).
2. **Conduction** — heat spreads through the solid graphite from hot spots (fuel) to
   cooler spots. Graphite is a decent conductor; that's why it's the matrix.
3. **Convection** — helium flowing through the coolant channels carries heat away at
   the channel walls. This is the *forced cooling*.
4. **Radiation** — any hot surface glows and sheds heat as thermal radiation. Here it
   matters at the **outer surface** of the block, which radiates to the reactor
   vessel. Radiation grows extremely fast with temperature (it scales with
   temperature to the 4th power), so it becomes the dominant escape route when things
   get hot.

The temperature settles where **heat made = heat removed** at every point. Hottest
place is inside the fuel; coolest is near the coolant channels and the outer edge.

**Why the 1600 °C number is the whole point:** the fuel is TRISO — tiny fuel kernels
each wrapped in ceramic shells that trap radioactive fission products. Those shells
stay intact up to about **1600 °C**. Keep the fuel below that and the reactor is safe
by design. So "peak fuel temperature vs 1600 °C" is *the* number this whole model
exists to produce.

> **Say it like this:** "The safety case for this fuel is a temperature limit —
> 1600 °C. My model's job is to predict the peak fuel temperature and show how much
> margin there is to that limit, in normal operation and in an accident."

---

## 3. The two scenarios (this is your headline story)

**Normal operation.** Full fission power, helium actively cooling the channels.
Result: peak fuel ≈ **822 °C**, which is **778 °C of margin** to the limit. Healthy.

**Loss of forced cooling (the interesting one).** Imagine the coolant blower fails.
Two things change:
- The reactor shuts down, so the fuel now only makes **decay heat** — leftover heat
  from radioactive decay, roughly a few percent of full power.
- With no flow, the coolant channels stop removing heat. The *only* way out is
  conduction to the outer surface, which then **radiates** the heat to the vessel.

Result: the block settles at ≈ **360 °C** — *cooler* than normal, because decay heat
is small and radiation easily handles it. Margin to the limit is now **1240 °C**.

> **Why cooler, not hotter? (be ready for this)** "It's counterintuitive, but the
> power dropped by about 30x the moment it shut down — it's only decay heat now. A
> small heat load radiating off a large surface reaches a low equilibrium
> temperature. That's exactly the passive-safety argument: you can walk away and it
> won't overheat."

This *is* the microreactor selling point, and it's a great thing to be able to
explain.

---

## 4. What "modeling it with FEM" actually means (no derivations)

You don't need to derive anything. Here's the honest conceptual picture:

- We can't solve the temperature everywhere with a single formula because the shape
  is complicated (holes, multiple materials). So we **chop the cross-section into
  thousands of little triangles** (the *mesh*). Corners of triangles are *nodes*.
- The **finite element method** is a systematic recipe for turning "heat is balanced
  everywhere" into "**a big system of equations for the temperature at every
  node**." Solve the system → temperature at every node → the whole field.
- **Boundary conditions** are how we tell the edges to behave: coolant walls remove
  heat by convection; the outer surface loses heat by convection + radiation.

That's it. FEM = mesh the shape, build a big equation system, solve for nodal
temperatures.

**One wrinkle — why it needs to iterate.** Two things make this problem *nonlinear*
(meaning you can't solve it in one shot):
1. Graphite's conductivity **changes with temperature**.
2. Radiation depends on **temperature to the 4th power**.

So the solver **guesses a temperature field, measures how far off the heat balance
is (the "residual"), corrects, and repeats** until the imbalance is essentially zero.
That repeat-until-balanced loop is **Newton's method**. The `convergence.png` figure
is literally that error shrinking to near-zero in about 4 steps.

> **Say it like this:** "It's nonlinear because conductivity varies with temperature
> and radiation goes as T-to-the-fourth, so the solver iterates — guess, check the
> imbalance, correct — until it converges. It takes about four iterations."

If they push into the *math* of how the correction is computed (the "Jacobian" /
"tangent"): see the honesty box in the interview prep doc. Short version — you
understand *what* it does (it's what makes the iteration converge fast); the detailed
derivation was AI-assisted and you can point them to the code and docs.

---

## 5. How we know the numbers are trustworthy (verification)

This is a strength you *can* own, because it's about logic, not derivation:

- **Test on a problem with a known answer.** We ran the solver on a made-up problem
  whose exact answer we know (the "Method of Manufactured Solutions"). The error
  shrank at exactly the rate the theory predicts as we refined the mesh (a measured
  "order of accuracy" of **1.999**, where 2.0 is the textbook target). That's strong
  evidence the code solves the equations correctly.
- **Energy audit.** We check that total heat generated equals total heat leaving.
  It balances to **0.000%**. If the code had a bug, energy wouldn't balance.
- **Two independent methods.** The plan is to reproduce a piece of this in OpenFOAM
  (a different, industry-standard method) and confirm they agree.

> **Say it like this:** "I didn't just trust the output — I verified it. It reproduces
> known-answer problems at the theoretically expected accuracy, energy balances to
> essentially zero, and I'm cross-checking against OpenFOAM. That's the same
> verification logic a safety analysis uses."

---

## 6. The files, and what each one is for

Run everything from inside the `fem_thermal/` folder.

| File | What it is | You'd say... |
|---|---|---|
| `materials.py` | The physical inputs: how conductivity varies with temperature, the two scenarios, the 1600 °C limit. | "where I set material properties and the two load cases" |
| `mesh.py` | Builds the hexagonal block geometry and chops it into triangles. | "the geometry and mesh" |
| `fem.py` | The solver — sets up and solves the temperature field, iterating until converged. | "the heat-transfer solver" |
| `verify.py` | The trust checks (known-answer tests, energy balance). | "my verification suite" |
| `main.py` | Runs both scenarios and makes the figures + results table. | "the driver that produces the results" |
| `plot_mesh.py` | Draws the mesh so you can see the geometry. | "a mesh viewer" |

**To reproduce the results live:**
```bash
cd fem_thermal
python verify.py     # prints the trust checks (should say 5/5 passed)
python main.py       # runs both scenarios, writes the figures
python plot_mesh.py  # draws the mesh
```

**The figures you'd show:**
- `figures/mesh.png` — the geometry (fuel = orange, coolant = blue, outer = red).
- `figures/temperature_fields.png` — the two temperature maps with peak-fuel marked.
- `figures/convergence.png` — the solver error shrinking (proof it converged).

---

## 7. Glossary (quick reference)

- **TRISO fuel** — fuel as tiny kernels wrapped in protective ceramic layers; safe up
  to ~1600 °C. The whole safety case rests on this number.
- **Prismatic block** — a solid graphite block with drilled channels for fuel and
  coolant; a common gas-cooled reactor fuel form.
- **Decay heat** — residual heat from radioactive decay after shutdown; a few percent
  of full power, and it can't be turned off. Passive cooling must handle it.
- **Conduction / convection / radiation** — the three heat-transfer modes: through a
  solid / carried by a moving fluid / emitted by a hot surface.
- **Boundary condition (BC)** — a rule for how heat behaves at an edge (e.g. "coolant
  removes heat here").
- **Mesh / node / element** — the shape chopped into small pieces (elements = the
  triangles, nodes = their corners).
- **FEM (finite element method)** — the recipe that turns the heat-balance physics
  into a solvable system of equations on the mesh.
- **Nonlinear / Newton iteration** — because properties depend on temperature, the
  solver must guess-and-correct repeatedly until balanced.
- **Residual** — the leftover heat imbalance; the solver drives it to ~zero.
- **Verification** — proving the code solves the equations correctly (vs reality,
  which is *validation*).
- **Margin** — how far the predicted peak temperature sits below the 1600 °C limit.
- **Peak fuel temperature** — the single most important output.

---

## 8. What to be honest about

You built this to **learn** gas-cooled reactor thermal analysis, and the
implementation was **AI-assisted**. That's a fine and increasingly normal thing to
say. What makes it land well:

- You can explain the **physics**, the **scenarios**, the **results**, and **why the
  verification matters** — all in this guide, all true.
- You don't claim to have hand-derived the finite element math. If asked, you say so
  plainly and redirect to what you *do* understand.

Interviewers are evaluating whether you can *think* about thermal problems and *use*
the tools — not whether you memorized a Jacobian derivation. Owning the story
honestly beats bluffing every time.

Next: read `07_interview_prep.md` for the likely questions and how to answer them.
