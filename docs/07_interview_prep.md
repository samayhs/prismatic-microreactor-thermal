# Interview Prep — Questions & How to Answer Them

**Document ID:** IP-001
Pair this with the study guide (SG-001). Answers are written the way *you* would say
them. Practice them out loud. The goal is calm, honest, specific.

---

## The 30-second pitch (lead with this if asked "tell me about a project")

> "I built a thermal model of a helium-cooled microreactor fuel block — the class of
> reactor Radiant builds. It predicts the peak fuel temperature, which is the key
> safety number because the TRISO fuel has to stay under about 1600 °C. I ran two
> cases: normal operation, and a loss-of-cooling accident where the block sheds heat
> only by radiation. The headline result is that it stays passively safe — about
> 1240 °C of margin even with no active cooling. And I verified the model against
> known-answer problems so the numbers are trustworthy, not just plausible."

---

## Tier 1 — Questions you should nail (physics, story, results)

**Q: What did you model, and why?**
> The temperature field in a gas-cooled reactor fuel block. The point was to find the
> peak fuel temperature and its margin to the 1600 °C TRISO limit, in normal
> operation and in a loss-of-cooling accident.

**Q: Why does the 1600 °C limit matter?**
> TRISO fuel is kernels wrapped in ceramic layers that hold in fission products up to
> ~1600 °C. Staying below that is the fundamental safety case, so peak fuel
> temperature is the number that matters most.

**Q: Walk me through the two scenarios.**
> Normal operation is full power with helium cooling the channels — peak fuel comes
> out around 822 °C. The accident case is loss of forced cooling: the reactor scrams
> so it's only making decay heat, and with no coolant flow the heat leaves only by
> radiating off the outer surface. It settles around 360 °C.

**Q: Wait — the accident case is *cooler*? That seems wrong.**
> It surprised me too at first, but it's correct. At shutdown the power drops roughly
> 30x — it's only decay heat. A small heat load radiating off a big surface reaches a
> low equilibrium temperature. That's the whole passive-safety idea: nothing has to
> work, and it still doesn't overheat.

**Q: What are the boundary conditions?**
> Coolant channel walls remove heat by convection — I model that with a heat-transfer
> coefficient. The outer surface loses heat by convection plus radiation to the
> vessel. In the accident case I turn the coolant off so radiation carries the load.

**Q: What did you conclude?**
> The design has large thermal margin in both cases, and it's passively safe — the
> accident case is well under the limit with no active systems. That matches the
> value proposition of a microreactor like Kaleidos.

---

## Tier 2 — Questions you can handle with the study guide

**Q: Why finite element method here?**
> The geometry is irregular — a block with drilled fuel and coolant channels and two
> materials. FEM handles complex geometry and mixed boundary conditions well by
> breaking the domain into a mesh and solving for temperature at each node.

**Q: What's the difference between FEM and CFD / finite volume?**
> They're both ways to turn the physics into a solvable system on a mesh. Finite
> element (what my Python solver uses) and finite volume (what OpenFOAM uses) discretize
> differently — finite volume enforces conservation cell-by-cell, which is why it's the
> standard for fluid flow. I'm using FEM for the solid conduction and OpenFOAM for the
> coolant flow, and cross-checking them against each other.

**Q: Why is the problem nonlinear?**
> Two reasons: graphite's thermal conductivity changes with temperature, and radiation
> scales with temperature to the fourth power. So the solver can't do it in one shot —
> it iterates: guess the field, measure the heat imbalance, correct, repeat until it
> converges. About four iterations here.

**Q: How do you know your results are correct?**
> Three ways. I ran the solver on a manufactured problem with a known exact answer and
> the error shrank at the theoretically expected rate as I refined the mesh — a measured
> order of accuracy of 1.999 versus a target of 2. I check that total heat generated
> equals total heat removed, and it balances to essentially zero. And I'm reproducing
> part of it in OpenFOAM to confirm two independent methods agree.

**Q: What are the limitations?**
> It's a 2D cross-section, so I'm not capturing axial coolant heat-up. The material
> properties are representative rather than a qualified vendor dataset, so absolute
> numbers are indicative — the emphasis is on the method and the margin behavior. And
> it's steady-state; a real accident is a transient, which is a natural next step.

---

## Tier 3 — Questions beyond "run tools" (handle honestly, don't bluff)

These go into the math you didn't derive yourself. The move is: **answer the concept,
be honest about the depth, redirect to what you understand.** This reads as mature,
not weak.

**Q: Derive the weak form / explain the finite element assembly.**
> "I understand it conceptually — FEM multiplies the governing equation by test
> functions and integrates to turn the PDE into a system of equations for the nodal
> temperatures, with the boundary conditions coming in through the boundary integrals.
> I'll be honest: I leaned on references and AI assistance for the detailed derivation
> and assembly code rather than doing it from memory. What I own is the physics setup,
> the boundary conditions, the scenarios, and verifying the results. Happy to walk you
> through the code."

**Q: Explain the Jacobian / Newton tangent.**
> "At a concept level, it's what makes the iteration correct itself efficiently — the
> sensitivity of the heat imbalance to a change in temperature, including the
> temperature-dependent conductivity and the radiation term. I verified it's right by
> comparing it against a finite-difference version and by confirming the solver
> converges quadratically. Deriving it by hand is past where I'm fluent today — it was
> AI-assisted — but I understand what role it plays."

**Q: Did you write all this yourself?**
> "I designed the problem, set up the physics and scenarios, ran and interpreted
> everything, and verified it. I used AI assistance for a lot of the numerical
> implementation — I treated it like a very capable pair-programmer. I built this to
> learn gas-cooled reactor thermal analysis, and I'm being upfront that the deep FEM
> math isn't something I'd derive from scratch yet. It's on my learning path."

**The universal redirect** when a question goes past your depth:
> "That's past where I'm fluent, so I won't fake it — but here's what I do understand
> about it, and here's the part of this project I can defend in detail…" then steer to
> physics / results / verification.

---

## What they're really evaluating (keep this in mind)

For a thermal-modeling role, especially at a fast-moving hardware company, they mostly
want to know:
1. Do you understand heat transfer and reason about it physically? — **Yes, show it.**
2. Can you set up, run, and interpret models and trust-check results? — **Yes, this
   project is exactly that.**
3. Are you honest about what you know and eager to learn? — **Yes, say so.**

They are *not* mainly testing whether you can derive a finite element formulation on a
whiteboard. Don't let a Tier-3 question rattle you — nailing Tier 1 and 2 and handling
Tier 3 with honesty is a strong interview.

---

## Two good questions to ask *them* (shows you think like a thermal engineer)

- "For Kaleidos, how much thermal margin do you design to on peak TRISO temperature
  under the bounding loss-of-cooling case?"
- "Is decay-heat rejection in your design conduction-dominated to the vessel, or does
  radiation carry most of it — and how does that change during transport?"

---

## Night-before checklist

- [ ] Re-read the study guide's "Say it like this" boxes.
- [ ] Say the 30-second pitch out loud 3 times.
- [ ] Be able to open the three figures and narrate each in one sentence.
- [ ] Re-read Tier 3 — get comfortable with the honest redirect so it feels natural.
- [ ] Have `python verify.py` and `python main.py` ready to run live if asked.
