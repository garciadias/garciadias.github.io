<script setup>
import RevealDeck from '@/components/RevealDeck.vue';

// Paper figures live in /public so they are served at the site root.
const fig = (n) => `${import.meta.env.BASE_URL}presentations/fla3/x${n}.png`
// Base URL for non-figure assets (e.g. the threat-model meme gif).
const baseUrl = import.meta.env.BASE_URL
</script>

<template>
  <RevealDeck :options="{ center: true }">
    <!-- 1 · Title -->
    <section class="title-slide center">
      <div class="eyebrow">Paper walkthrough · <a href="https://arxiv.org/abs/2603.10063">arXiv 2603.10063 (submitted to IEEE Journal of Biomedical and Health Informatics on March 2026)</a></div>
      <h1>What can FLIP learn from FLAAA?</h1>
      <a href="https://arxiv.org/abs/2603.10063">
        <p class="subtitle">
          <strong>FLA<sup>3</sup></strong> — Federated Learning with <em>AAA</em> <br>
          (Authentication · Authorisation · Accounting)
        </p>
        <p class="byline">
          Fan Zhang et al. 2026<br />
          Cambridge / KCL / QMUL /  Maastricht University / Flower Labs  <br />
        </p>
      </a>
      <a href="https://garciadias.github.io/#/presentations/fla3-governance-federated-learning" class="link">
        <p class="venue-note">AMIGO team meeting · 2 June 2026</p><br>
        <p class="venue-note">Slides: tinyurl.com/amigo-fla3</p>
      </a>
      <aside class="notes">
        Frame the talk for our own team: this is a read of the FLA³ paper through one question —
        what can FLIP take from it? Their thesis, and the one to hold onto: in regulated healthcare
        the binding constraint isn't data leakage, it's unauthorised computation — and runtime
        governance is feasible at zero accuracy cost. The cyan FLIP asides throughout track where we
        already match them and where we don't.
      </aside>
    </section>

    <!-- 2 · The setup · abstract & introduction -->
    <section>
      <div class="eyebrow">The problem</div>
      <h2>The data exists — but it can't be pooled</h2>
      <div class="cols" style="--n: 3">
        <div class="panel">
          <h3>The need</h3>
          <p class="small">
            Clinical AI needs <strong>large, diverse, representative</strong> data — which lives
            scattered across institutions and borders.
          </p>
        </div>
        <div class="panel">
          <h3>The barrier</h3>
          <p class="small">
            <strong>GDPR · UK GDPR · HIPAA</strong> restrict or prohibit centralised aggregation
            of health data — especially across jurisdictions.
          </p>
        </div>
        <div class="panel fla3">
          <h3>The accepted answer</h3>
          <p class="small">
            <strong>Federated learning</strong> trains a shared model while raw data
            <em>stays local</em>. The data never moves.
          </p>
        </div>
      </div>
      <p>
        So the standard story ends here. But in <strong>regulated</strong> healthcare it doesn't:
        most FL frameworks remain <strong>proof-of-concept</strong>, assuming trusted participants
        and ideal conditions — and leaving the governance questions of
        <em>who may compute, on what, and for how long</em> unanswered.
      </p>
      <div class="flip-note">
        <span class="flip-tag">FLIP</span> Our own platform lives in exactly this world — FL across
        NHS Trusts. So as we go, the question to hold is: <em>where do FLA<sup>3</sup> and FLIP make
        the same call, and where do they diverge? What can we learn from them?</em>
      </div>
      <aside class="notes">
        The opener — set the scene, don't argue yet. Everyone knows FL keeps data local; the hook is
        the turn at the end: that's necessary but not sufficient once a regulator is in the room —
        who may compute, on what, for how long. That gap is the whole paper. Introduce FLIP here as
        the peer system we read against the rest of the way.
      </aside>
    </section>

    <!-- 3 · The gap (evidence) -->
    <section>
      <div class="eyebrow">Why this matters now</div>
      <h2>Governance is the systematic blind spot</h2>
      <div class="kpis" style="--n: 3">
        <div class="panel center">
          <div class="stat">~5%</div>
          <div class="stat-label">of healthcare FL studies reach real-world deployment</div>
        </div>
        <div class="panel center">
          <div class="stat">98%</div>
          <div class="stat-label">of 89 reviewed FL methods lacked node authentication</div>
        </div>
        <div class="panel center">
          <div class="stat">0</div>
          <div class="stat-label">offered a peer-reviewed, openly specified governance impl.</div>
        </div>
      </div>
      <p class="muted small">
        Differential privacy &amp; secure aggregation stop <em>information leakage</em> — they
        do not enforce <strong>institutional accountability</strong> or time-bounded
        authorisation. Those are orchestration-layer problems.
      </p>
      <aside class="notes">
        Let the numbers do the work: ~5% of healthcare FL reaches deployment, 98% lacked node
        authentication, and zero shipped an openly specified governance implementation. The line to
        voice: privacy tech (DP, secure aggregation) stops information leakage but not unauthorised
        computation — governance is an orchestration-layer problem most FL simply skips.
      </aside>
    </section>

    <!-- 4 · Existing frameworks fall short (+ FLIP aside) -->
    <section>
      <div class="eyebrow">Where current FL frameworks stop</div>
      <h2>Security ≠ governance</h2>
      <table>
        <thead>
          <tr><th>Framework</th><th>Provides</th><th>Missing at the orchestration layer</th></tr>
        </thead>
        <tbody>
          <tr>
            <td><strong>Flower</strong></td>
            <td>mTLS comms, basic user auth, lifecycle mgmt</td>
            <td>No study-scoped authz or time-based validity</td>
          </tr>
          <tr>
            <td><strong>PySyft</strong></td>
            <td>Secure aggregation, differential privacy</td>
            <td>Cryptographic privacy, not workflow governance</td>
          </tr>
          <tr>
            <td><strong>NVIDIA FLARE</strong></td>
            <td>Role-based authz, site access control</td>
            <td>Configured up front — no per-round/per-study runtime enforcement, no auto-expiry</td>
          </tr>
        </tbody>
      </table>
      <div class="flip-note">
        <span class="flip-tag">FLIP</span> Builds on <strong>NVFLARE + Flower</strong>. 
        FLIP adds project scoping, approval governance &amp; user auth at the hub and per-Trust isolation, 
        at its <em>hub/Trust APIs</em>, but not a per-round policy decision point inside the FL
        runtime, or automation on study-level policies.
      </div>
      <aside class="notes">
        Honest framing for an FL crowd — Flower, PySyft and FLARE are all good systems; the gap is
        runtime, per-study, time-bounded enforcement bound into orchestration. The note for us: FLIP
        builds on NVFLARE + Flower, so it inherits exactly this gap at the FL layer — we add our
        governance at the hub/Trust APIs instead.
      </aside>
    </section>

    <!-- 5 · FLA3 in one line (+ FLIP aside) -->
    <section>
      <div class="eyebrow">The proposal</div>
      <h2>A governance <span class="accent">control plane</span> inside FL</h2>
      <p>
        FLA<sup>3</sup> integrates <strong>XACML-compliant ABAC</strong>, cryptographic
        accounting and <strong>study-scoped federation</strong> directly into the orchestration
        layer — framework-agnostic, demonstrated by extending <strong>Flower</strong>.
      </p>
      <div class="cols" style="--n: 3">
        <div class="panel">
          <h3>Authentication</h3>
          <p class="small">mTLS, institutional certificates — identity only, never privilege.</p>
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov total">✅</span>
            Cognito <em>user</em> auth + <strong>institutional</strong> mTLS identity.
          </div>
        </div>
        <div class="panel">
          <h3>Authorisation</h3>
          <p class="small">XACML PDP — study-scoped, role-aware, time-bounded, fail-closed.</p>
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            Per-<em>project</em> Cognito roles at the hub — not re-checked
            <strong>per round</strong>, no temporal auto-expiry.
          </div>
        </div>
        <div class="panel">
          <h3>Accounting</h3>
          <p class="small">Tamper-evident JWS audit of every decision.</p>
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            Platform / API logs only — not <strong>cryptographically tamper-evident</strong>.
          </div>
        </div>
      </div>
      <div class="flip-note">
        <span class="flip-tag">FLIP</span> FLIP's design choices are mostly aligned with FLA<sup>3</sup>. However, 
        we could implement a per-round PDP and cryptographic audit in FLIP if we wanted. This can create unnecessary
        friction, but we definitely should revisit the trade-offs and consider adding these features in the future.
      </div>
      <aside class="notes">
        The one-line pitch: governance as a first-class control plane, not paperwork bolted on after.
        Walk the three AAA blocks — the FLIP verdict under each is the point: we match on
        authentication (✅), but only partially on authorisation (no per-round PDP) and accounting
        (no signed audit). Plant the idea that both are ours to close.
      </aside>
    </section>

    <!-- 6 · Five contributions (roadmap) -->
    <section>
      <div class="eyebrow">Contributions · what this paper delivers</div>
      <h2>Five contributions</h2>
      <ol class="contribs small">
        <li>
          <strong>Regulatory-derived governance requirements</strong> — a systematic
          cross-jurisdictional analysis (HIPAA · GDPR · DPDPA · ECOWAS) formalised as enforceable
          system properties <span class="req">R1</span>–<span class="req">R5</span>.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            We have not done a formal cross-jurisdictional analysis, we try to follow best practices and common principles, but we haven't mapped them to enforceable properties in the same way. 
          </div>
        </li>
        <li>
          <strong>Policy-driven AAA framework</strong> — XACML authorisation + cryptographic
          accounting on Flower's auth; <strong>fail-closed</strong> at every FL lifecycle point,
          with cryptographically signed audit records.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            Auth + project authorisation at the hub APIs — no XACML PDP or signed audit records, but we do follow fail-closed lifecycle checks.
          </div>
        </li>
        <li>
          <strong>Multi-study federation</strong> — each study is an independent unit with its own
          participant set, policies &amp; temporal validity; many concurrent studies on one
          platform.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov total">✅</span>
            Close match — FLIP already isolates many concurrent <em>projects</em> on one platform.
          </div>
        </li>
        <li>
          <strong>Governance-preserving personalisation</strong> — with <strong>FedMAP</strong> on
          the INTERVAL 25-centre data, federation beats per-site training and matches centralised
          — enforcement costs nothing.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov total">✅</span>
            Runs FL training (NVFLARE / Flower), we don't support FedMAP yet, but this is an aggregation strategy. We support arbitrary aggregation strategies, so FedMAP could be added without changing the governance layer.
          </div>
        </li>
        <li>
          <strong>Open-source reference implementation</strong> —
          <code>github.com/bloodcounts/FLAAA</code>, a reference architecture for governance-aware
          FL in regulated healthcare.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov total">🏆</span>
            Also open-source (<code>londonaicentre/FLIP</code>) — production code rather than a
            reference architecture.
          </div>
        </li>
      </ol>
      <aside class="notes">
        The paper's five contributions, each with our honest FLIP verdict. We already match on
        multi-study (✅) and open-source (🏆); the partials — formal cross-jurisdictional
        requirements, XACML PDP + signed audit, and FedMAP — are the candidate borrows. This doubles
        as the roadmap: 1–3 are the mechanism, 4 is the proof, 5 is what we could clone.
      </aside>
    </section>

    <!-- 7 · Five requirements -->
    <section>
      <div class="eyebrow">Contribution 1 · derived from HIPAA · GDPR · DPDPA · ECOWAS · ICH-GCP</div>
      <h2>Five enforceable requirements</h2>
      <ul class="reqs">
        <li>
          <span class="req">R1</span><strong>Authenticated institutional participation</strong> — identity bound to ethics/DSA attributes (approvals issued to <em>institutions</em>).
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            Users authenticated (Cognito) + per-Trust isolation — but identity isn't bound to ethics/DSA attributes as policy inputs.
          </div>
        </li>
        <li>
          <span class="req">R2</span><strong>Study-scoped authorisation</strong> — approval for one protocol never transfers to another (purpose limitation).
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov total">✅</span>
            Close match — every action is scoped to an approved project; authorisation never transfers across projects.
          </div>
        </li>
        <li>
          <span class="req">R3</span><strong>Role-based, least-privilege access</strong> — participant vs observer; no federation-level escalation.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov total">✅</span>
            Cognito roles give per-project least-privilege — no participant/observer split at the round level.
          </div>
        </li>
        <li>
          <span class="req">R4</span><strong>Temporal validity</strong> — authorisation auto-expires when approvals/DSAs lapse.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            The real gap — approvals are checked at onboarding, not auto-expired per round when a DSA lapses.
          </div>
        </li>
        <li>
          <span class="req">R5</span><strong>Accounting &amp; auditability</strong> — every security-relevant action recorded for attribution.
          <div class="flip-inline">
            <span class="flip-tag">FLIP</span><span class="cov partial">🟡</span>
            Actions logged at the platform / API level — but not cryptographically signed or tamper-evident.
          </div>
        </li>
      </ul>
      <aside class="notes">
        R1–R5 recur across every jurisdiction they analysed — enforceable system properties, not
        box-ticking checked once at onboarding. Per requirement, flag where FLIP stands: R2 and R3
        are close matches; R4 (temporal auto-expiry) is the real gap; R1 and R5 are partial. These
        are the same five we return to in the discussion.
      </aside>
    </section>

    <!-- 8 · Threat model -->
    <section>
      <div class="eyebrow">Threat model</div>
      <h2>Authenticated ≠ authorised</h2>
      <div class="fig-split" style="--cols: 1.7fr 1fr; align-items: center">
        <div>
          <p class="small muted">
            Adversary = valid credentials, network access &amp; ability to invoke federation APIs —
            but may lack valid clinical approval at request time. Four governance-layer threats:
          </p>
          <div class="cols" style="--n: 2">
            <div class="panel"><h3><span class="req">A1</span> Unauthorised participation</h3><p class="small">No/expired approval, or out of protocol scope → <span class="req">R1</span><span class="req">R2</span><span class="req">R4</span></p></div>
            <div class="panel"><h3><span class="req">A2</span> Privilege misuse</h3><p class="small"> Authorised user goes beyond their role → <span class="req">R3</span><span class="req">R2</span></p></div>
            <div class="panel"><h3><span class="req">A3</span> Accountability evasion</h3><p class="small">Change logs to avoid responsibility → <span class="req">R5</span></p></div>
            <div class="panel"><h3><span class="req">A4</span> Policy bypass</h3><p class="small">Omit attributes (incomplete context) to exploit default-permit → <strong class="pink">fail-closed</strong></p></div>
          </div>
        </div>
        <div class="panel center">
          <img :src="`${baseUrl}presentations/fla3/security-meme.gif`" alt="Mr. Robot main character looking suspicious." style="min-height: 500px;" />
        </div>
      </div>
      <aside class="notes">
        The mental-model flip: stop picturing the hoodie-in-a-server-room hacker (the meme) — the
        real adversary is authenticated, holds valid credentials and API access, and simply isn't
        authorised at that moment. A1–A4 are all governance-layer. Flag the scope boundary now:
        Byzantine / malicious-update defence is out of scope and complementary — it comes back on
        the limitations slide.
      </aside>
    </section>

    <!-- 9 · Architecture (Figure 1) -->
    <section>
      <div class="eyebrow">System architecture · Fig. 1</div>
      <h2>Three layers, governance in the middle</h2>
      <div class="fig-split" style="--cols: 1fr 0.85fr">
        <div class="figure">
          <img :src="fig(1)" alt="FLA3 system architecture — SuperLink, SuperNodes, ClientApps with AAA governance layer" style="max-height: 460px" />
        </div>
        <div>
          <ul class="small">
            <li><strong>SuperLink</strong> — central coordination; per-study ServerApps, aggregation rounds.</li>
            <li><strong>SuperNode</strong> — site-local gateway at each hospital.</li>
            <li><strong>ClientApp</strong> — ephemeral, study-scoped train/eval processes.</li>
            <li>AAA box wraps the server: <strong>PDP</strong> + cryptographic audit.</li>
          </ul>
          <div class="flip-note">
            <span class="flip-tag">FLIP</span> Maps to FLIP's <em>Central Hub API</em> +
            per-Trust <em>Trust/imaging/data-access APIs</em> over XNAT/Orthanc/OMOP — but FLIP's
            control points are the APIs, not a PDP woven into FL rounds.
          </div>
        </div>
      </div>
      <aside class="notes">
        Three layers, governance in the middle: SuperLink coordinates, SuperNodes sit at each site,
        ClientApps are ephemeral per-study processes, and the AAA box (PDP + audit) wraps the server.
        For us: this maps onto FLIP's Central Hub + per-Trust APIs — but our control points are the
        APIs, not a PDP woven into the FL rounds.
      </aside>
    </section>

    <!-- 10 · Multi-study federation + deployment -->
    <section>
      <div class="eyebrow">Multi-study federation &amp; deployment</div>
      <h2>Many studies, one platform, no inbound ports</h2>
      <div class="fig-split" style="--cols: 1.7fr 1fr; align-items: center">
        <div>
          <ul class="small">
            <li><strong>Client-initiated gRPC unary</strong> calls — no inbound connectivity needed; deployable behind restrictive hospital firewalls &amp; TREs (outbound-only).</li>
            <li><strong>Concurrent studies</strong> — each with its own member set, policies &amp; validity window; identified by <strong>immutable study IDs</strong>.</li>
            <li>Credential compromise is <strong>contained to a single study</strong>.</li>
          </ul>
          <div class="flip-note">
            <span class="flip-tag">FLIP</span> FLIP solves the same network limitations with a similar client-initiated design.
            FLIP's per-Trust APIs are also scoped to projects, which provide the same isolation and guarantees as FLA3's study-scoped design.
          </div>
        </div>
        <div class="panel center">
          <img :src="`${baseUrl}presentations/fla3/flex-meme.gif`" alt="Kid flexing, captioned FLEX" style="max-height: 320px;" />
        </div>
      </div>
      <aside class="notes">
        Two real-world wins worth dwelling on: client-initiated gRPC needs no inbound ports, so it
        deploys behind hospital firewalls and TREs; and study-scoping contains a credential
        compromise to a single study. The flex is earned — and FLIP earns it too: our per-Trust,
        project-scoped, client-initiated design gives the same guarantees.
      </aside>
    </section>

    <!-- 11 · Authentication -->
    <section>
      <div class="eyebrow">AAA · Authentication</div>
      <h2>Identity, not privilege</h2>
      <p>
        Mutual TLS with <strong>institutional credentials</strong> — every node presents a
        certificate issued to a recognised healthcare institution before connecting. Integrated
        with Flower authentication; certificate attributes feed the PDP.
      </p>
      <p class="muted">Authentication binds <em>who you are</em>. It grants <strong>nothing</strong>.</p>
      <div class="flip-note">
        <span class="flip-tag">FLIP</span> FLIP authenticates users via <strong>AWS
        Cognito</strong> and institutions with TLS certificates. This provides the same level of identity assurance as FLA<sup>3</sup>,
        but also extends the identity to individual users, which allows for more granular access control and auditing.
      </div>
      <aside class="notes">
        The key distinction to land: authentication is identity only — it grants nothing. mTLS
        institutional certs just feed attributes to the PDP; authorisation happens next. For us this
        is a ✅, and arguably finer-grained: FLIP authenticates institutions with TLS and individual
        users via Cognito, which is what gives us per-user audit.
      </aside>
    </section>

    <!-- 12 · Authorisation: the request model (core) -->
    <section>
      <div class="eyebrow">AAA · Authorisation — the core mechanism</div>
      <h2>One request, three predicates, fail-closed</h2>
      <p class="small muted">
        PDP = <strong>Luas</strong>, an XACML-compliant engine. Flower server components act as
        Policy Enforcement Points. 
        <br> Request = ⟨ node, action (train | eval), study, timestamp ⟩.
        <br> ⚠️ The <strong>Luas</strong> package
        actual implementation is not open-source. 
      </p>
      <ul class="checklist small">
        <li><strong>Study validity</strong> — now &lt; study expiry</li>
        <li><strong>Membership validity</strong> — node is a current, unexpired member of <em>that</em> study</li>
        <li><strong>Role eligibility</strong> — train ⇒ <em>participant</em>; eval ⇒ <em>participant</em> or <em>observer</em></li>
      </ul>
      <p class="small">
        <strong class="green">PERMIT</strong> only if all three hold. Enforced at <em>node
        activation</em>, <em>study join</em>, and <strong>every aggregation round</strong>.
        Missing attribute or eval error → <strong class="pink">DENY</strong>.
      </p>
      <aside class="notes">
        The core mechanism: one request ⟨node, action, study, time⟩, three predicates, PERMIT only if
        all hold — and it's re-checked every aggregation round, which is what makes temporal validity
        real (an approval lapsing mid-training blocks the very next round, not at some future audit).
        Be honest about the caveat on the slide: the Luas PDP engine itself isn't open-source.
      </aside>
    </section>

    <!-- 13 · Obligation → mechanism mapping (Table 1) -->
    <section>
      <div class="eyebrow">Table 1 · Obligations → mechanisms</div>
      <h2>Regulation, mapped to enforcement</h2>
      <table>
        <thead>
          <tr><th>Obligation</th><th>Regulatory basis</th><th>FL mechanism</th><th>Component</th></tr>
        </thead>
        <tbody>
          <tr><td>Institutional identity</td><td>HIPAA §164.312; GDPR Art.5(2)</td><td>mTLS + institutional creds</td><td>Flower auth + PDP attrs <span class="req">R1</span></td></tr>
          <tr><td>Institutional approval</td><td>HRA/REC, EU CTR, ICMR</td><td>Approval/DSA status as attributes</td><td>PDP membership validation</td></tr>
          <tr><td>Purpose limitation</td><td>GDPR Art.5(1)(b); DPDPA</td><td>Study-scoped policies, immutable IDs</td><td>Study validation <span class="req">R2</span></td></tr>
          <tr><td>Access control</td><td>HIPAA §164.312; Caldicott</td><td>Role &amp; action enforcement</td><td>PDP role policies <span class="req">R3</span></td></tr>
          <tr><td>Temporal validity</td><td>Ethics approval validity</td><td>Time-bounded conditions, auto-expiry</td><td>PDP temporal predicates <span class="req">R4</span></td></tr>
          <tr><td>Auditability</td><td>HIPAA §164.312(b); GDPR Art.30</td><td>Integrity-protected audit logs</td><td>JWS audit logging <span class="req">R5</span></td></tr>
        </tbody>
      </table>
      <aside class="notes">
        The spine of the contribution: every legal obligation maps to a concrete, runtime-checkable
        mechanism and a named component (R1–R5). Don't read the table out — point at the
        through-line: regulation on the left, an enforceable mechanism on the right, for each row.
      </aside>
    </section>

    <!-- 14 · Accounting -->
    <section>
      <div class="eyebrow">AAA · Accounting</div>
      <h2>Tamper-evident by construction</h2>
      <ul class="small">
        <li>Every authorisation decision → a structured <strong>JSON</strong> audit record.</li>
        <li>Signed with a <strong>JSON Web Signature</strong> (<code>ES256</code>) — any edit invalidates the signature.</li>
        <li>Binds <em>subject · action · study ID · timestamp · policy version</em> → non-repudiation for audit &amp; forensics.</li>
      </ul>
      <div class="flip-note">
        <span class="flip-tag">FLIP</span> FLIP relies on platform/API logging. Cryptographic,
        signed, tamper-evident records (R5) are an FLA<sup>3</sup> strength worth importing for
        regulatory defensibility.
      </div>
      <aside class="notes">
        Tamper-evident by construction: every decision is a JWS-signed record, so you can't quietly
        rewrite history — the direct answer to the A3 accountability-evasion threat. For us this is
        the clearest thing to import: FLIP logs at the platform/API level, but those logs aren't
        cryptographically signed.
      </aside>
    </section>

    <!-- 15 · Validation: operational deployment -->
    <section>
      <div class="eyebrow">Validation · Operational feasibility</div>
      <h2>5 institutions · 4 countries · 4 legal regimes</h2>
      <p>
        Platform infrastructure deployed across <strong>BloodCounts! Consortium</strong> sites in
        the <strong>UK, Netherlands, India &amp; The Gambia</strong> — GDPR, DPDPA (consent-centric)
        and ECOWAS-aligned regimes; heterogeneous network postures (some outbound-only).
      </p>
      <p class="muted small">
        Governance enforcement operated correctly under realistic network &amp; regulatory
        constraints. <em>This study proves it deploys — the next proves it predicts.</em>
      </p>
      <aside class="notes">
        Be precise: this is the operational study — does governance actually hold up live across 4
        countries and 4 legal regimes (some outbound-only)? It proves the platform deploys. The
        clinical numbers next come from a separate, UK simulated federation — don't conflate the two.
      </aside>
    </section>

    <!-- 16 · Clinical setup: what is actually being measured -->
    <section>
      <div class="eyebrow">Validation · Clinical utility · INTERVAL trial (ISRCTN24760606)</div>
      <h2>The test: flag <span class="accent">iron deficiency</span> from a routine blood count</h2>
      <p class="small">
        Binary task — predict <strong>iron deficiency</strong> (serum ferritin &lt; 15&nbsp;µg/L) from
        <strong>routine full-blood-count parameters only</strong>, no extra blood draw. Scored by
        <strong>ROC–AUC</strong> (1.0 = perfect, 0.5 = chance). Identical model throughout: a small
        2-layer MLP (64 → 32). <strong>54,446</strong> samples · <strong>25</strong> UK donation
        centres of <em>very</em> uneven size (CV 64%) — realistic deployment.
      </p>
      <div class="cols" style="--n: 3">
        <div class="panel">
          <h3>Individual</h3>
          <p class="small">Each centre trains on <em>its own data only</em> — today's status quo.</p>
        </div>
        <div class="panel fla3">
          <h3>Federated · FedMAP</h3>
          <p class="small">Personalised FL across all 25 centres, <strong>governance enforced</strong>. What FLA<sup>3</sup> delivers.</p>
        </div>
        <div class="panel">
          <h3>Centralised</h3>
          <p class="small">Pool every centre's data in one place — an <em>upper bound</em> that governance makes infeasible in practice.</p>
        </div>
      </div>
      <p class="muted small" style="margin-top: 0.4em">
        The question for the next two slides: <strong>can governed federation close the gap to the
        (forbidden) centralised ceiling?</strong>
      </p>
      <aside class="notes">
        Spend time here — it frames both result slides. Three ways to train: individual (status quo),
        federated FedMAP (governance on), and centralised (the upper bound governance forbids). The
        whole question: does governed federation reach that forbidden ceiling? ROC–AUC ≈ the
        probability the model ranks a random iron-deficient donor above a random healthy one.
      </aside>
    </section>

    <!-- 17 · Clinical result: no accuracy cost (Fig 2) -->
    <section>
      <div class="eyebrow">Result 1 · Predictive performance · Fig. 2</div>
      <h2>Governance at <span class="green">zero accuracy cost</span></h2>
      <div class="fig-split" style="--cols: 1fr 0.95fr">
        <div class="figure">
          <img :src="fig(2)" alt="Per-centre ROC-AUC: individual training vs FedMAP, with centralised reference line" style="max-height: 420px" />
        </div>
        <div>
          <p class="small muted" style="margin-top: 0">
            Each point = one of 25 centres · box = spread across centres ·
            <strong class="pink">red dashed line</strong> = centralised upper bound.
          </p>
          <div class="kpis" style="--n: 1; gap: 0.4em">
            <div class="panel center"><div class="stat">0.872</div><div class="stat-label">FedMAP mean ROC–AUC — lands <em>exactly</em> on the centralised ceiling</div></div>
            <div class="panel center"><div class="stat">+0.027</div><div class="stat-label">vs individual training (0.845 → 0.872; 95% CI 0.019–0.036)</div></div>
          </div>
          <p class="small muted">
            Wilcoxon p &lt; 0.001 · Cohen's d = 1.28 — <strong>federation reaches the centralised
            ceiling without pooling any data.</strong>
          </p>
        </div>
      </div>
      <aside class="notes">
        Headline for an ML audience: you don't trade accuracy for governance. FedMAP's mean (0.872)
        lands exactly on the centralised dashed line, well above individual training (0.845) — and it
        gets there without pooling a single record. That's the whole "zero cost" claim in one chart.
      </aside>
    </section>

    <!-- 18 · Clinical result: equity (Fig 3 + Fig 4) -->
    <section>
      <div class="eyebrow">Result 2 · Equity across centres · Figs. 3–4</div>
      <h2>The weakest sites gain the most</h2>
      <div class="cols" style="--n: 2; align-items: start">
        <div class="figure center">
          <img :src="fig(3)" alt="Federation gain vs baseline performance, Pearson r = -0.74" style="max-height: 300px" />
          <p class="small muted" style="margin: 0.3em 0 0">
            <strong>Fig. 3</strong> — each centre's <em>gain</em> (FedMAP − individual) vs its
            starting score. Downward trend (<strong>r = −0.74</strong>): the worse you started, the
            more you gain.
          </p>
        </div>
        <div class="figure center">
          <img :src="fig(4)" alt="Regional median ROC-AUC maps: individual vs FedMAP" style="max-height: 300px" />
          <p class="small muted" style="margin: 0.3em 0 0">
            <strong>Fig. 4</strong> — regional median ROC–AUC, same colour scale (individual → FedMAP).
            Federation lifts the map and <strong>evens it out</strong>.
          </p>
        </div>
      </div>
      <p class="small" style="margin-top: 0.3em">
        <strong>24/25 (96%)</strong> centres improved · inter-centre SD <strong>0.029 → 0.020</strong>
        · regional spread <strong>0.117 → 0.065</strong> (≈ −44%).
      </p>
      <aside class="notes">
        The fairness story: federation doesn't just lift the mean (Fig. 2), it compresses the spread.
        The worse a centre started, the more it gained (r = −0.74), and regional variation nearly
        halved. The one centre that dipped already had the best baseline — least to gain. This is the
        equity argument regulators care about for multi-centre AI.
      </aside>
    </section>

    <!-- 18b · The worry scale (interlude) -->
    <section>
      <div class="eyebrow">Are we being replaced by FLA<sup>3</sup>?</div>
      <h2>So… how worried should we be?</h2>
      <div class="cols" style="--n: 3; align-items: end">
        <div class="panel center">
          <img :src="`${baseUrl}presentations/fla3/worry-1.gif`" alt="Calm dog basking in the sun" style="min-height: 230px; border-radius: 8px; margin: 0.3em 0" />
        </div>
        <div class="panel center">
          <img :src="`${baseUrl}presentations/fla3/worry-3.gif`" alt="Fry squinting suspiciously" style="min-height: 330px; border-radius: 8px; margin: 0.3em 0" />
        </div>
        <div class="panel center">
          <img :src="`${baseUrl}presentations/fla3/worry-5.gif`" alt="Person hyperventilating into a paper bag" style="min-height: 330px; border-radius: 8px; margin: 0.3em 0" />
        </div>
      </div>
      <div class="scale-bar"></div>
      <div class="scale-ticks"><span>1</span><span>2</span><span>3</span><span>4</span><span>5</span></div>
      <aside class="notes">
        Comic beat with a real question — "are we (FLIP) being replaced by FLA³?" Answer over the
        next slides: no. Read the scale left to right: 1 = calm (data stays local and it works),
        3 = vigilant (authenticated ≠ authorised), 5 = panic (an approval lapsed mid-study and
        nothing stopped the next round). The point under the joke: runtime governance slides you back
        toward 1. Sets up the honest limitations.
      </aside>
    </section>

    <!-- 19 · Limitations (honest) -->
    <section>
      <div class="eyebrow">Honest limitations</div>
      <h2><span class="pink">Limitations</span></h2>
      <ul class="crosslist small">
        <li>Clinical evaluation is a <strong>simulated</strong> federation — not yet live end-to-end training across the deployed sites.</li>
        <li>No risk profile for <strong>Byzantine</strong> participants with valid credentials submitting malicious updates.</li>
        <li>No risk profile for arbitrary code execution.</li>
        <li>No consideration of <strong>patient opt-out status</strong>.</li>
        <li>Composition with <strong>differential privacy / secure aggregation</strong> not yet empirically validated.</li>
        <li>XACML policies are <strong>manually authored</strong> — room for specification error.</li>
      </ul>
      <div class="flip-note">
        <span class="flip-tag">FLIP</span> Several of these are exactly where <strong>we're</strong>
        strong, but we also have limitations to overcome.
      </div>
      <aside class="notes">
        Honest about what the paper doesn't yet claim: simulated (not live) federation, no risk
        profile for Byzantine participants or arbitrary code execution, no patient opt-out handling,
        and manually authored policies. Voicing these builds trust — and several are exactly where
        FLIP is strong (we run live, multi-Trust). The gaps run both ways; that's the bridge to the
        discussion.
      </aside>
    </section>

    <!-- 20 · Discussion / thanks -->
    <section>
      <div class="eyebrow">Discussion · AMIGO team meeting · 2 June 2026</div>
      <h2>What should <span class="cyan">FLIP</span> take from FLA<sup>3</sup>?</h2>
      <div class="cols" style="--n: 4">
        <div class="panel flip">
          <h3>Governance Requirements</h3>
          <p class="small">They have done a important cross-jurisdictional healthcare governance analysis we can certainly learn from.</p>
        </div>
        <div class="panel flip">
          <h3>Per-round PDP?</h3>
          <p class="small">Re-check authorisation every aggregation round — worth the latency at NHS scale, or is hub-level enough for our threat model?</p>
        </div>
        <div class="panel flip">
          <h3>Temporal auto-expiry?</h3>
          <p class="small">Tie participation to DSA / ethics validity so it lapses automatically — not only checked at onboarding.</p>
        </div>
        <div class="panel flip">
          <h3>Signed audit?</h3>
          <p class="small">Cryptographically tamper-evident decision records, for regulatory defensibility.</p>
        </div>
      </div>
      <p class="center accent" style="margin-top: 0.7em">Thank you — let's discuss.</p>
      <p class="center muted small" style="margin-top: 0.2em">
        Paper: <span class="accent">arXiv 2603.10063</span> · FLA<sup>3</sup>:
        <span class="accent">github.com/bloodcounts/FLAAA</span> · FLIP:
        <span class="cyan">github.com/londonaicentre/FLIP</span>
      </p>
      <aside class="notes">
        Close on our turf — four things FLA³ suggests for FLIP, framed as open questions, not
        prescriptions: their cross-jurisdictional governance analysis to learn from, a per-round PDP,
        temporal auto-expiry, and signed audit. Open the floor on per-round PDP vs latency first —
        that's the live debate for us at NHS scale.
      </aside>
    </section>

    <!-- 21 · Appendix: acronym reference -->
    <section>
      <div class="eyebrow">Appendix · for reference</div>
      <h2>Acronyms</h2>
      <div class="glossary">
        <div class="gl-group">
          <h4>Concepts &amp; systems</h4>
          <p><b>FL</b> — Federated Learning</p>
          <p><b>AAA</b> — Authentication, Authorisation &amp; Accounting</p>
          <p><b>FLA³</b> — <a href="https://github.com/bloodcounts/FLAAA" class="link">Federated Learning with AAA</a></p>
          <p><b>FLIP</b> — <a href="https://github.com/londonaicentre/FLIP" class="link">FL and Interoperability Platform</a></p>
          <p><b>FedMAP</b> — Federated Maximum A Posteriori</p>
          <p><b>MLP</b> — Multi-Layer Perceptron</p>
        </div>
        <div class="gl-group">
          <h4>Access control &amp; crypto</h4>
          <p><b>ABAC</b> — Attribute-Based Access Control</p>
          <p><b>XACML</b> — eXtensible Access Control Markup Language</p>
          <p><b>PDP</b> — Policy Decision Point</p>
          <p><b>PEP</b> — Policy Enforcement Point</p>
          <p><b>TLS / mTLS</b> — (mutual) Transport Layer Security</p>
          <p><b>JWS</b> — JSON Web Signature</p>
          <p><b>JSON</b> — JavaScript Object Notation</p>
          <p><b>ES256</b> — ECDSA P-256 + SHA-256 (JWS algorithm)</p>
        </div>
        <div class="gl-group">
          <h4>Platforms, infra &amp; standards</h4>
          <p><b>API</b> — Application Programming Interface</p>
          <p><b>gRPC</b> — gRPC Remote Procedure Call</p>
          <p><b>AWS</b> — Amazon Web Services</p>
          <p><b>NVFLARE</b> — NVIDIA FL Application Runtime Environment</p>
          <p><b>TRE</b> — Trusted Research Environment</p>
          <p><b>XNAT</b> — eXtensible Neuroimaging Archive Toolkit</p>
          <p><b>OMOP</b> — Observational Medical Outcomes Partnership (data model)</p>
          <p><b>DICOM</b> — Digital Imaging &amp; Communications in Medicine</p>
        </div>
        <div class="gl-group">
          <h4>Regulation &amp; governance</h4>
          <p><b>GDPR / UK GDPR</b> — (UK) General Data Protection Regulation</p>
          <p><b>HIPAA</b> — Health Insurance Portability &amp; Accountability Act</p>
          <p><b>DPDPA</b> — Digital Personal Data Protection Act (India)</p>
          <p><b>ECOWAS</b> — Economic Community of West African States</p>
          <p><b>DSA</b> — Data Sharing Agreement</p>
          <p><b>HRA / REC</b> — Health Research Authority / Research Ethics Committee</p>
          <p><b>EU CTR</b> — EU Clinical Trials Regulation</p>
          <p><b>ICMR</b> — Indian Council of Medical Research</p>
          <p><b>ICH-GCP</b> — Int'l Council for Harmonisation – Good Clinical Practice</p>
          <p><b>NHS</b> — National Health Service</p>
        </div>
        <div class="gl-group">
          <h4>Clinical &amp; statistics</h4>
          <p><b>FBC</b> — Full Blood Count</p>
          <p><b>ROC–AUC</b> — Receiver Operating Characteristic – Area Under the Curve</p>
          <p><b>CV</b> — Coefficient of Variation</p>
          <p><b>CI</b> — Confidence Interval</p>
          <p><b>DP</b> — Differential Privacy</p>
          <p><b>ISRCTN</b> — Int'l Standard Randomised Controlled Trial Number</p>
        </div>
      </div>
      <aside class="notes">
        Reference only — don't walk through it. Leave it up during Q&amp;A so any acronym is one
        glance away.
      </aside>
    </section>
  </RevealDeck>
</template>
