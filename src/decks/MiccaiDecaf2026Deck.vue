<script setup>
import RevealDeck from '@/components/RevealDeck.vue';
</script>

<template>
  <RevealDeck :options="{ center: true }">
    <!-- 1 · Title -->
    <section class="title-slide center">
      <div class="eyebrow">Working draft · MICCAI 2026 · DECAF workshop</div>
      <h1>Federated chest X-ray learning, UK ⇄ Thailand</h1>
      <p class="subtitle">
        Cross-continental federated fine-tuning of a chest X-ray
        <strong>foundation model</strong> — a proof of concept built on <span class="cyan">FLIP</span>.
      </p>
      <p class="byline">
        Anonymised submission · FLIP team<br />
        UK academic institutions ⇄ private hospital group, Thailand
      </p>
      <p class="venue-note">Internal draft preview · 2026</p>
    </section>

    <!-- 2 · Introduction -->
    <section>
      <div class="eyebrow">Introduction</div>
      <h2>Diverse data, locked in silos</h2>
      <div class="cols" style="--n: 3">
        <div class="panel">
          <h3>The need</h3>
          <p class="small">
            Imaging AI improves with <strong>large, diverse, multi-institutional</strong> data —
            spread across hospitals and borders which comply with <strong>privacy, governance &amp; security</strong> 
            rules across jurisdictions.
          </p>
        </div>
        <div class="panel">
          <h3>The barrier</h3>
          <p class="small">
            Setting a FL network up is <strong>not trivial</strong>, taking months of engineering per project — and
            almost never reported across the <em>academic ⇄ private-sector</em> boundary. Most FL research focuses on
            the incremental algorithmic improvements, not on building robust federated networks that can be reused
            across projects.
          </p>
        </div>
        <div class="panel fla3">
          <h3>Paper contribution</h3>
          <p class="small">
            This is the first report of a <strong>cross-continental, cross-sector</strong> FL deployment linking UK
            academic institutions with a private hospital group in Thailand — and the first to characterise the
            <strong>engineering cost</strong> of standing up a federation that can be reused for future projects.
          </p>
        </div>
      </div>
      <p>
        We show here a prof of concept for federated fine-tuning of a chest X-ray foundation model, with all data local and only model weights exchanged.
      </p>
      <aside class="notes">
        Set the scene: everyone accepts FL keeps data local. The hook is the turn — the binding cost
        isn't the algorithm, it's the repeated engineering of connecting heterogeneous institutions,
        and almost nobody reports a federation that crosses the academic/private-sector line. That gap
        is what this paper attacks.
      </aside>
    </section>

    <!-- 3 · The idea -->
    <section>
      <div class="eyebrow">The idea</div>
      <h2>Not a new algorithm — reusable <span class="cyan">infrastructure</span></h2>
      <div class="cols" style="--n: 2">
        <div class="panel flip">
          <h3>Reuse FLIP</h3>
          <p class="small">
            FLIP — the Federated Learning &amp; Interoperability Platform — is built to cut per-project
            engineering overhead, and is <strong style="color: red;">already validated?</strong> at <strong>population scale across NHS
            hospital networks</strong>.
          </p>
          <p class="muted small">
            This is a <strong>deployment study and implementation roadmap</strong> — feasibility,
            reproducibility and initial performance — not a new federated-learning method. The work showcase the platform
            and invite the cummunity to contribute and reuse it for future federated medical-AI projects.
          </p>

        </div>
        <div class="panel">
          <h3>Push it somewhere new</h3>
          <p class="small">
            We demonstrate FLIP's portability by deploying it across a <strong>cross-continental, cross-institution</strong>
            federation — UK academic institutions ⇄ a private hospital group in Thailand — and characterise the
            <strong>governance, technical &amp; operational barriers</strong> encountered along the way and contrast it
            with the traditional centralised approach. We also demonstrate the results on a controlled federated fine-tuning 
            experiment of a chest X-ray foundation model, comparing zero-shot, local-only, federated and centralised
            training under matched conditions.
          </p>
        </div>
      </div>
      <aside class="notes">
        The framing the whole paper rests on: we are not proposing a new FL algorithm. We are testing
        whether existing infrastructure (FLIP), proven inside the NHS, transfers to a harder setting —
        across continents and across the academic/private-sector boundary. Keep expectations honest:
        scope is feasibility and a roadmap, not clinical deployment.
      </aside>
    </section>

    <!-- 4 · Research question & hypothesis -->
    <section>
      <div class="eyebrow">Research question &amp; hypothesis</div>
      <h2>Can FLIP carry the load?</h2>
      <div class="cols" style="--n: 2">
        <div class="panel">
          <h3>Question</h3>
          <p class="small">
            Can FLIP <strong>simplify deployment</strong> of a federated network across large,
            heterogeneous <em>and cross-sector</em> infrastructures while delivering measurable
            research value in a cross-continental imaging task?
          </p>
        </div>
        <div class="panel fla3">
          <h3>Hypothesis</h3>
          <p class="small">
            A UK ⇄ Thailand FLIP workflow trains a shared chest X-ray model with <strong>all data
            local</strong>; the federated model <strong>beats local-only</strong> and
            <strong>approaches centralised</strong> under matched conditions.
          </p>
          <p class="small">
            The <strong>operational cost</strong> of standing up the federation is lower than
            building a bespoke workflow from scratch, and the <strong>engineering burden</strong>
            is characterised for future FLIP roll-outs.
          </p>
          <p class="small">
            FLIP meet all security, governance and privacy requirements for cross-institutional collaboration,
            reducing the friction of IT approvals, firewall configuration, credential management, cloud setup and monitoring.
          </p>

        </div>

      </div>
      <p class="muted small">
        As a deployment study, we additionally characterise the <strong>institutional and technical
        barriers</strong> encountered — to inform future FLIP roll-outs.
      </p>
      <aside class="notes">
        Two halves: the engineering question (does FLIP lower the connection cost across a cross-sector,
        cross-continental boundary?) and the modelling hypothesis (federated &gt; local-only,
        federated ≈ centralised). The barrier-characterisation is a deliberate third deliverable — the
        failures and friction are findings, not footnotes.
      </aside>
    </section>

    <!-- 5 · Contributions -->
    <section>
      <div class="eyebrow">Contributions · what this paper delivers</div>
      <h2>Four contributions</h2>
      <ol class="contribs small">
        <li>
          <strong>First deployment of FLIP beyond the NHS</strong> — a live cross-continental,
          cross-sector federation linking UK academic institutions with a private hospital group in
          Thailand.
        </li>
        <li>
          <strong>Institutional &amp; technical barriers identified</strong> from the deployment, with
          concrete implications for future cross-sector FLIP roll-outs.
        </li>
        <li>
          <strong>A reproducible federated fine-tuning workflow</strong> for a pretrained chest X-ray
          model, with a controlled comparison of zero-shot, local-only, federated &amp; centralised
          training under matched conditions.
        </li>
      </ol>
      <aside class="notes">
        Four contributions, and notice none of them is "a better aggregator." 1 is the deployment
        first; 2 is the controlled four-arm comparison; 3 is the insistence on system metrics next to
        AUC; 4 turns the deployment pain into a reusable roadmap. This slide doubles as the outline of
        the rest of the talk.
      </aside>
    </section>

    <!-- 6 · Dive into -->
    <section>
      <div class="eyebrow">Deeper dive</div>
      <h2>What we build on</h2>
      <div class="cols" style="--n: 3">
        <div class="panel">
          <h3>FL in medical imaging</h3>
          <p class="small">
            Federated performance hinges on <strong>client heterogeneity</strong> — sample size,
            prevalence, protocol, device, labels. Aggregation often <em>narrows inter-site variance</em>
            even when it doesn't beat a single-site baseline.
          </p>
        </div>
        <div class="panel">
          <h3>Deployment &amp; comms</h3>
          <p class="small">
            Real systems must handle bandwidth, latency, update size, synchronisation, aggregation time
            and dropped connections — <strong>acute cross-continentally</strong>, where big updates cost
            time and money.
          </p>
        </div>
        <div class="panel">
          <h3>Data heterogeneity</h3>
          <p class="small">
            Non-Identically distributed (non-IID) data is the norm in medical imaging, and cross-continental
            federation is a <strong>distribution-shift</strong> problem. Thai-vs-Western show case ethinical, protocol
            and device heterogeneity.
          </p>
        </div>
      </div>
      <p class="muted small">
        Governance adds a per-project burden of its own — which is why most prominent FL applications are
        <strong>single-project</strong> rather than reusable, multi-project platforms.
      </p>
      <aside class="notes">
        Three threads. First: heterogeneity is the dominant variable in FL imaging, and aggregation's
        real win is often variance reduction, not raw accuracy. Second: communication is a first-class
        constraint here because the link is intercontinental. Third: heterogeneity is plural — we care
        specifically about Thai-vs-Western distribution shift. The governance line motivates why a
        platform (FLIP) beats one-off pipelines.
      </aside>
    </section>

    <!-- 7 · FLIP architecture -->
    <section>
      <div class="eyebrow">FLIP architecture &amp; workflow</div>
      <h2>One hub, two nodes, weights only on the wire</h2>
      <div class="cols" style="--n: 3; align-items: stretch">
        <div class="panel">
          <h3>UK node</h3>
          <p class="small">
            Starter package · <strong>OMOP</strong> (EHR) · <strong>XNAT</strong> (images) · Flare
            client. Western synthetic cohort.
          </p>
        </div>
        <div class="panel fla3 center">
          <h3>Central Hub · AWS</h3>
          <p class="small">
            Federated-learning server (<strong>Flare</strong>) · periodic
            <strong>FedAvg or FedProx</strong> weight averaging · privacy-preserving monitoring.
          </p>
        </div>
        <div class="panel">
          <h3>Thailand node</h3>
          <p class="small">
            Same starter package · OMOP · XNAT · Flare client. Thai synthetic cohort.
          </p>
        </div>
      </div>
      <p class="center small muted" style="margin-top: 0.5em">
        ⟵ model weights only (FedAvg/FedProx) ⟶ &nbsp;·&nbsp; raw images <strong>never leave a node</strong>
      </p>
      <p class="small">
        Per-site nodes install from a <strong>starter package</strong> to minimise engineering overhead;
        FL backend is NVIDIA Flare (Flower is also supported). Each site ingests data via SQL, loads the
        <em>same</em> pretrained model + fine-tuning code, trains locally, then averages across sites.
      </p>
      <aside class="notes">
        Walk it left-to-right: identical nodes at each site (OMOP for structured data, XNAT for images,
        a Flare client), an AWS hub running the Flare server and FedAvg. The one sentence that matters
        — only weights cross the wire, images never do. The starter package is the whole point: it's how
        FLIP collapses the per-project setup cost that Related Work flagged.
      </aside>
    </section>

    <!-- 8 · Study design -->
    <section>
      <div class="eyebrow">Methods · study design</div>
      <h2>A technical proof of concept</h2>
      <ul class="checklist small">
        <li><strong>Two nodes</strong> represent the UK and Thailand sites.</li>
        <li>An <strong>AWS</strong> environment provides coordination &amp; model aggregation.</li>
        <li>Raw images <strong>remain at their node</strong> throughout federated training.</li>
        <li>Only <strong>model weights / updates</strong> are transmitted.</li>
        <li>Designed as a <strong>feasibility study</strong>, not a breakthrough clinical trial — results are illustrative.</li>
      </ul>
      <aside class="notes">
        Keep this short — it's the contract for everything that follows. Two sites, AWS-mediated
        aggregation, images pinned to their node, weights-only on the wire, and an explicit "this is a
        PoC, not a clinical trial" disclaimer so nobody over-reads the results.
      </aside>
    </section>

    <!-- 9 · Model -->
    <section>
      <div class="eyebrow">Methods · model</div>
      <h2>Fine-tune, don't train from scratch</h2>
      <ul class="small">
        <li>Start from a <strong>pretrained chest X-ray foundation model / backbone</strong>.</li>
        <li><strong>Ark+</strong> is the preferred model — subject to confirming code availability, input
          requirements &amp; FLIP-pipeline compatibility.</li>
        <li>Adapt the pretrained representation to the target <strong>binary classification</strong> task —
          never trained from scratch.</li>
        <li>Fine-tuning strategy (frozen / partially / fully trainable encoder) <strong>fixed before</strong>
          the primary experiments.</li>
        <li><strong>Same initialisation &amp; comparable optimisation budgets</strong> across all arms — for a
          fair comparison.</li>
      </ul>
      <aside class="notes">
        Two things to land: (1) we adapt a foundation model rather than training from scratch — that's
        what makes a 2-site PoC viable at all; (2) fairness is engineered in — same init, same budget,
        strategy frozen up front — so any difference between arms is the federation, not the
        hyperparameters. Ark+ is the front-runner but still provisional.
      </aside>
    </section>

    <!-- 10 · Data: synthetic by design -->
    <section>
      <div class="eyebrow">Methods · data</div>
      <h2>No real radiographs train the model</h2>
      <div class="cols" style="--n: 2">
        <div class="panel">
          <h3>UK client · Western</h3>
          <p class="small">
            Frontal CXRs generated <strong>text-to-image</strong> (GPT-image), prompts specifying each
            target finding. A western-style synthetic cohort, conditioned on <em>no real image</em>.
          </p>
        </div>
        <div class="panel">
          <h3>Thailand client · Thai</h3>
          <p class="small">
            Touches just <strong>5 real radiographs</strong> (one prototypical case per finding) →
            ≈ 1,000 synthetic via image editing:
          </p>
          <p class="small" style="margin: 0.25em 0 0">
            <strong>① Exemplar-conditioned editing</strong> — diversify variants that preserve the
            class-defining pathology.<br />
            <strong>② Normalise-then-inject</strong> — edit to a no-finding image, then inject the target
            pathology for controlled paired normal/abnormal cases.
          </p>
        </div>
      </div>
      <p class="muted small">
        The 5 Thai exemplars are <strong>held out from, and disjoint with, the evaluation set</strong> —
        no real test image or its derivatives appear in training. Minimal real-data footprint by design.
      </p>
      <aside class="notes">
        This is the unusual part — neither client trains on real images. The UK side is pure
        text-to-image; the Thai side leans on exactly five purposively chosen real exemplars expanded two
        ways: edit-to-diversify, and normalise-then-inject for controlled paired cases. The privacy
        argument is the held-out-disjoint guarantee: those five never touch evaluation, directly or via
        their derivatives.
      </aside>
    </section>

    <!-- 11 · Diagnostic-alignment filtering -->
    <section>
      <div class="eyebrow">Methods · data quality control</div>
      <h2>A VLM judge gates every synthetic image</h2>
      <p class="small">
        Generators don't guarantee the finding they were prompted for — so before any image enters
        training, an <strong>automated alignment gate</strong> runs <em>locally</em> at each client
        (exchanging no data).
      </p>
      <div class="cols" style="--n: 2">
        <div class="panel">
          <h3>The judge</h3>
          <p class="small">
            A vision-language model (<strong>Gemini</strong>) decides: a plausible frontal CXR? the
            intended finding present &amp; consistent with the label? Only <strong>confirmed-aligned</strong>
            images are kept — clients <em>over-generate</em> to fill the pool. Each verdict carries a short
            <strong>justification + localisation</strong> for auditability.
          </p>
        </div>
        <div class="panel">
          <h3>Validating the judge</h3>
          <p class="small">
            A general VLM isn't a specialist reader, and generator/judge errors may correlate. So we
            validate against a <strong>clinician</strong> on a stratified sample (percent agreement,
            <strong>Cohen's κ</strong>) and report <strong>per-class rejection rates</strong> — watching
            pneumothorax &amp; rib fracture, where both models are likeliest to err.
          </p>
        </div>
      </div>
      <aside class="notes">
        Synthetic data is only useful if it actually shows the labelled finding, so a VLM judge filters
        every image locally — accept only if it's a plausible CXR and the target finding is present.
        Crucially we don't trust the judge blindly: we benchmark it against a human reader with Cohen's
        kappa and watch the hard findings (pneumothorax, rib fracture) where generator and judge errors
        could correlate.
      </aside>
    </section>

    <!-- 12 · Partitioning & heterogeneity -->
    <section>
      <div class="eyebrow">Methods · partitioning &amp; heterogeneity</div>
      <h2>Western at one node, Thai at the other</h2>
      <ul class="small">
        <li>Primary federated config: <strong>synthetic Western → UK node</strong>,
          <strong>synthetic Thai → Thailand node</strong>.</li>
        <li><strong>Centralised reference</strong> combines both datasets in one place — for benchmarking only.</li>
        <li><strong>Local-only baselines</strong> trained independently at each node on its own data.</li>
        <li>Train / validation / test splits <strong>non-overlapping &amp; fixed</strong> across comparable experiments.</li>
        <li>Report class balance, sample count, prevalence &amp; image characteristics per partition.</li>
        <li>Evaluate both <strong>in-distribution</strong> and <strong>cross-distribution</strong> (Thai ⇄ Western) generalisation.</li>
      </ul>
      <aside class="notes">
        The partition is the experiment: Western data on one node, Thai on the other, so the federation
        genuinely spans a distribution shift. Centralised pools both (benchmark only); local-only is each
        site alone. The split is frozen across arms, and we test both in-distribution and cross-distribution
        — the latter is where federation should earn its keep.
      </aside>
    </section>

    <!-- 13 · Experimental comparisons -->
    <section>
      <div class="eyebrow">Methods · experimental comparisons</div>
      <h2>Four models, matched conditions</h2>
      <div class="cols" style="--n: 2">
        <div class="panel">
          <h3>Zero-shot · 𝓜<sub>zs</sub></h3>
          <p class="small">Pretrained model evaluated on Thai &amp; Western test sets — <em>no</em> fine-tuning.</p>
        </div>
        <div class="panel">
          <h3>Local-only · 𝓜<sub>local</sub></h3>
          <p class="small">Fine-tuned on a <em>single</em> site; each local model evaluated on <strong>both</strong> sites' tests.</p>
        </div>
        <div class="panel fla3">
          <h3>Federated · 𝓜<sub>FL</sub></h3>
          <p class="small">Trained across both sites via <strong>FedAvg/FedProx</strong>, alternating local &amp; global rounds, monitored from the hub.</p>
        </div>
        <div class="panel">
          <h3>Centralised · 𝓜<sub>ref</sub></h3>
          <p class="small">Trained on <strong>joined</strong> data — the reference. Nothing but the training set differs from 𝓜<sub>FL</sub>.</p>
        </div>
      </div>
      <p class="muted small">
        𝓜<sub>FL</sub> vs 𝓜<sub>ref</sub> is the clean test: identical everything <em>except</em> whether
        data was pooled.
      </p>
      <aside class="notes">
        The four-arm design is the heart of the evaluation. Zero-shot sets the floor; local-only shows
        what a site gets alone (and how it generalises to the other); centralised is the upper-bound
        reference. The key contrast is federated vs centralised — they differ in one thing only, whether
        the data was pooled — so the gap between them is the price (or not) of keeping data local.
      </aside>
    </section>

    <!-- 14 · Evaluation -->
    <section>
      <div class="eyebrow">Evaluation</div>
      <h2>Accuracy <em>and</em> the cost of getting there</h2>
      <div class="cols" style="--n: 2">
        <div class="panel">
          <h3>Predictive performance</h3>
          <ul class="small" style="margin-top: 0.2em">
            <li><strong>Primary:</strong> ROC–AUC.</li>
            <li><strong>Secondary:</strong> accuracy, sensitivity, specificity, F1, PR-AUC.</li>
            <li>Reported per distribution (Thai / Western) <em>and</em> pooled.</li>
            <li><strong>Client-level</strong> results — does FL narrow or widen site gaps?</li>
          </ul>
        </div>
        <div class="panel">
          <h3>System &amp; deployment</h3>
          <ul class="small" style="margin-top: 0.2em">
            <li>Communication rounds · per-round local time · aggregation time · total wall-clock.</li>
            <li>Update size per client/round · total communication volume.</li>
            <li>Connection failures, failed rounds, retries &amp; recovery time.</li>
            <li>Hardware utilisation · <strong>manual engineering effort</strong> to install &amp; run.</li>
          </ul>
        </div>
      </div>
      <aside class="notes">
        The contribution-3 slide in practice: we refuse to report AUC alone. The right column is what
        most FL papers omit — rounds, wall-clock, update size, total bytes, failures and recovery, and
        crucially the human engineering effort. For a deployment study, those are the result, as much as
        the accuracy is.
      </aside>
    </section>

    <!-- 15 · Reproducibility -->
    <section>
      <div class="eyebrow">Evaluation · reproducibility &amp; quality control</div>
      <h2>Built to be re-run</h2>
      <ul class="checklist small">
        <li>Fixed <strong>random seeds</strong> for partitioning, initialisation &amp; training.</li>
        <li>Key experiments <strong>repeated across seeds</strong> where feasible — report mean &amp; variability, not the best run.</li>
        <li>Log software versions, hardware, FLIP config, instance type, optimiser, LR, batch size, local epochs, rounds &amp; stopping rules.</li>
        <li>Verify all arms use <strong>identical test sets</strong> &amp; equivalent preprocessing.</li>
        <li>A structured experiment log links every result to its config &amp; checkpoint.</li>
      </ul>
      <aside class="notes">
        Reproducibility is a stated deliverable, so it gets its own slide. The honest detail is "report
        mean and variability, not the best run" — and that every reported number is traceable to a config
        and checkpoint. This is what lets the next site re-run the study rather than reinvent it.
      </aside>
    </section>

    <!-- 16 · Results placeholder -->
    <section>
      <div class="eyebrow">Results · draft</div>
      <h2>What we will report</h2>
      <div class="cols" style="--n: 4">
        <div class="panel">
          <h3>Governance</h3>
            <p class="small muted">Ethical considerations, regulatory compliance, data governance, and stakeholder engagement.</p>
        </div>
        <div class="panel">
          <h3>Model performance</h3>
          <p class="small muted">Zero-shot, local-only, centralised &amp; federated AUC — per distribution, client-level &amp; cross-distribution, with uncertainty.</p>
        </div>
        <div class="panel">
          <h3>FL training behaviour</h3>
          <p class="small muted">Training/validation curves over rounds, convergence, between-client variability, stability across seeds.</p>
        </div>
        <div class="panel">
          <h3>System &amp; deployment</h3>
          <p class="small muted">Install/config effort, communication volume, latency, aggregation &amp; wall-clock time, failures, recovery &amp; lessons.</p>
        </div>
      </div>
      <p class="center muted small" style="margin-top: 0.6em">
        <em>Numbers to be completed after experimentation &amp; deployment — this is a working draft.</em>
      </p>
      <aside class="notes">
        Be upfront that this is a draft: the results are scaffolded, not filled. Use this slide to show
        the team exactly which tables and figures are coming, so the discussion can be about whether we're
        measuring the right things — not about missing numbers.
      </aside>
    </section>

    <!-- 17 · Discussion -->
    <section>
      <div class="eyebrow">Discussion · expected</div>
      <h2>What we hope to show</h2>
      <div class="cols" style="--n: 2">
        <div class="panel">
          <h3>Did the federation work?</h3>
          <p class="small">
            Whether FLIP supported <strong>end-to-end cross-continental FL with no raw-data transfer</strong>,
            and whether 𝓜<sub>FL</sub> beat local-only, approached centralised, and
            <strong>narrowed cross-site gaps</strong>.
          </p>
        </div>
        <div class="panel flip">
          <h3>Did FLIP lower the cost?</h3>
          <p class="small">
            Whether FLIP <strong>reduced the engineering burden</strong> of standing up the workflow —
            approvals, firewalls, credentials, cloud config, storage &amp; monitoring — versus building it
            from scratch.
          </p>
        </div>
      </div>
      <p class="muted small">
        Distinguish <strong>predictive gains</strong> from the <strong>operational value</strong> of a
        repeatable international workflow — both are outcomes.
      </p>
      <aside class="notes">
        Two questions map to the two halves of the hypothesis. Left: the modelling result — did federation
        help, and did it help the weaker side most? Right: the engineering result — did FLIP actually cut
        the setup cost? Insist they're separable wins: even flat accuracy with a dramatically cheaper
        deployment is a positive result for this paper.
      </aside>
    </section>

    <!-- 18 · Limitations -->
    <section>
      <div class="eyebrow">Honest limitations</div>
      <h2><span class="pink">Limitations</span></h2>
      <ul class="crosslist small">
        <li>A <strong>two-site</strong> proof of concept — may not represent larger federations.</li>
        <li><strong>Synthetic</strong> images don't reproduce all properties of real clinical data; artefacts &amp; biases may shape results.</li>
        <li>Primary task limited to <strong>binary</strong> classification.</li>
        <li>Single <strong>foundation model</strong>.</li>
        <li>The model is performing quite strongly even without fine-tuning.</li>
        <li>Successful deployment ≠ clinical effectiveness, safety, fairness or regulatory readiness.</li>
        <li>No model deployment to real-world settings.</li>
      </ul>
      <aside class="notes">
        State the limits before anyone else does — it buys credibility. The big ones: two sites, fully
        synthetic data, a single binary task, one model and one aggregator, and the disclaimer that a
        working deployment says nothing yet about clinical safety. These also seed the future-work slide.
      </aside>
    </section>

    <!-- 19 · Future work
    <section>
      <div class="eyebrow">Future work</div>
      <h2>Where this goes next</h2>
      <ul class="checklist small">
        <li>Extend to <strong>more institutions</strong> &amp; more heterogeneous hardware.</li>
        <li>Evaluate <strong>private clinical datasets</strong> after governance, security &amp; ethics approvals.</li>
        <li>Compare more <strong>aggregation, personalisation &amp; communication-efficient</strong> strategies.</li>
        <li>Move to <strong>multi-label</strong> chest X-ray tasks &amp; clinically relevant external validation.</li>
        <li>Build automated <strong>deployment checks, monitoring, failure recovery &amp; reproducibility tooling</strong> into FLIP.</li>
      </ul>
      <aside class="notes">
        Future work falls into two buckets: scale (more sites, real data, multi-label, external
        validation) and hardening FLIP itself (automated deployment checks, monitoring, recovery). The
        second bucket is where this PoC's barrier-findings feed straight back into the platform.
      </aside>
    </section> -->

    <!-- 20 · Conclusion
    <section class="center">
      <div class="eyebrow">Conclusion</div>
      <h2>A feasibility study, not a new algorithm</h2>
      <p>
        We evaluate whether <span class="cyan">FLIP</span> can support <strong>reproducible
        cross-continental federated fine-tuning</strong> of a chest X-ray foundation model between the UK
        and Thailand — weighing <strong>model performance</strong> against the
        <strong>operational realities</strong> of federated deployment.
      </p>
      <p class="muted small">
        The contribution is a practical feasibility study, implementation record &amp; roadmap for future
        international medical-AI collaboration.
      </p>
      <p class="accent" style="margin-top: 0.6em">Thank you — let's discuss.</p>
      <aside class="notes">
        Land the plane on the thesis: the value here is the deployment and the roadmap, not a novel
        method. Invite the team to push on two things in discussion — is the synthetic-data design
        defensible, and are we measuring the right deployment costs? Those are the live questions for the
        camera-ready.
      </aside>
    </section>
 -->
    <!-- 21 · Acronyms appendix -->
    <!-- <section>
      <div class="eyebrow">Appendix · for reference</div>
      <h2>Acronyms</h2>
      <div class="glossary">
        <div class="gl-group">
          <h4>Learning &amp; platform</h4>
          <p><b>FL</b> — Federated Learning</p>
          <p><b>FLIP</b> — <a href="https://github.com/londonaicentre/FLIP" class="link">FL &amp; Interoperability Platform</a></p>
          <p><b>FedAvg</b> — Federated Averaging</p>
          <p><b>PoC</b> — Proof of Concept</p>
          <p><b>Ark+</b> — chest X-ray foundation model</p>
          <p><b>VLM</b> — Vision-Language Model</p>
        </div>
        <div class="gl-group">
          <h4>Clinical &amp; imaging</h4>
          <p><b>CXR</b> — Chest X-ray</p>
          <p><b>EHR</b> — Electronic Health Record</p>
          <p><b>OMOP</b> — Observational Medical Outcomes Partnership (data model)</p>
          <p><b>XNAT</b> — eXtensible Neuroimaging Archive Toolkit</p>
          <p><b>NHS</b> — National Health Service</p>
        </div>
        <div class="gl-group">
          <h4>Metrics &amp; stats</h4>
          <p><b>ROC–AUC</b> — Receiver Operating Characteristic – Area Under the Curve</p>
          <p><b>PR-AUC</b> — Precision–Recall AUC</p>
          <p><b>F1</b> — F1 score</p>
          <p><b>κ</b> — Cohen's kappa (inter-rater agreement)</p>
          <p><b>LR</b> — Learning Rate</p>
        </div>
        <div class="gl-group">
          <h4>Infra &amp; governance</h4>
          <p><b>AWS</b> — Amazon Web Services</p>
          <p><b>gRPC</b> — gRPC Remote Procedure Call</p>
          <p><b>SQL</b> — Structured Query Language</p>
          <p><b>GDPR</b> — General Data Protection Regulation</p>
          <p><b>HIPAA</b> — Health Insurance Portability &amp; Accountability Act</p>
          <p><b>DSA</b> — Data Sharing Agreement</p>
        </div>
      </div>
      <aside class="notes">
        Reference only — leave it up during Q&amp;A so any acronym is one glance away.
      </aside>
    </section> -->
  </RevealDeck>
</template>
