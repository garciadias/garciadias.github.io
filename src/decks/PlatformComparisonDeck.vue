<script setup>
import RevealDeck from "@/components/RevealDeck.vue";
</script>

<template>
    <RevealDeck :options="{ center: true }">
        <!-- 1 · Title. Framed as a methods talk: the audience already knows FLIP,
         so the interesting content is how the comparison was made defensible. -->
        <section class="title-slide center">
            <div class="eyebrow">Federated Learning Platform Comparison</div>
            <h1>FLIP in its ecosystem</h1>
            <p class="subtitle">
                Inclusion criteria, a reproducible search, <br />
                and a 30-platform capability audit scored on
                <em>code, not claims</em>
            </p>
            <p class="byline">
                Rafael Garcia-Dias · AI Centre for Value-Based Healthcare<br />
                King's College London
            </p>
            <p class="venue-note">AMIGO team meeting · 15th September 2026</p>
            <aside class="notes">
                (~30 s) Frame it: this is the methods talk behind the table; how
                we decided who is in the comparison, how we screened them, and
                what we found. The punchline: each of the four searches has a
                different blind spot, and the only way to cover the full field
                is to use all four together and say precisely which one caught
                what.
            </aside>
        </section>

        <!-- 2 · Where this comes from. -->
        <section>
            <div class="eyebrow">They asked for a comparison</div>
            <h2>Why we needed a comprehensive platform comparison?</h2>
            <ul class="dotlist small">
                <li>
                    Position FLIP accurately in the field. No re-implementing or
                    re-inventing the wheel.
                </li>
                <li>
                    MICCAI workshop reviewers asked for a comprehensive
                    comparison for a workshop paper.
                </li>
            </ul>
            <div class="slide-body">
                <div class="cols" style="--n: 3">
                    <div class="panel flip">
                        <h3>Capability</h3>
                        <p class="small">
                            30 platforms &times; 9 criteria. What each one
                            <strong>ships</strong>, scored by reading a local
                            clone at a recorded commit. The paper only where no
                            source is reachable.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Project health</h3>
                        <p class="small">
                            Active contributors, commits, PRs, latest release.
                            Whether each one is
                            <strong>still alive</strong>. Metric families follow
                            CHAOSS.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Deployment</h3>
                        <p class="small">
                            What each has actually
                            <strong>demonstrated</strong>, sites, jurisdictions,
                            modality, scale.
                        </p>
                    </div>
                </div>
                <p class="small" style="margin-top: 0.7em">
                    Main deliverables: Inclusion Criteria well defined and 3
                    evidence base tables for our positioning claim.
                </p>
            </div>
            <aside class="notes">
                (~1 min) Quick orientation for anyone who hasn't read the draft.
                Three tables. The capability table is the contested one because
                it's what backs our novelty claim. Stress the scoring rule early
                Or rather, the current public source rather than the paper;
                because everything later depends on it.
            </aside>
        </section>

        <!-- 3 · The reviewer's question. The honest admission that motivates the
         whole talk. -->
        <section>
            <div class="eyebrow">The problem</div>
            <h2>How to select the platforms for comparison?</h2>
            <div class="slide-body">
                <h3>
                    We started with the ones we knew, but moved to a systematic
                    search to refine the set.
                </h3>
                <div class="cols" style="--n: 2">
                    <div class="panel">
                        <h3>Failure mode one</h3>
                        <p class="small">
                            We
                            <strong
                                >miss a platform that already does what we claim
                                is novel</strong
                            >.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Failure mode two</h3>
                        <p class="small">
                            We include only platforms that make us look good.
                            Nobody has to allege bad faith,
                            <strong
                                >an unstated selection rule is
                                unfalsifiable</strong
                            >, which is enough.
                        </p>
                    </div>
                </div>
                <div class="flip-note">
                    <span class="flip-tag">FLIP</span> The comparison is only
                    useful to us if we can be rigorous and critical of our own
                    position. We want to present a comparison that can be
                    reproduced, verified, and challenged.
                </div>
            </div>
            <aside class="notes">
                (~1 min 15) This is the slide where I admit the gap. Say it
                plainly: the selection was expert judgement, which is how nearly
                every comparison table in this field is built, and it's not good
                enough for a paper making a novelty claim. Both failure modes
                matter, but the second is the one that kills you in review,
                because you can't disprove it without a stated rule.
            </aside>
        </section>

        <!-- 4 · Unit of analysis. CODA is the case that forces the rule. -->
        <section>
            <div class="eyebrow">Criteria · 1 of 4</div>
            <h2>A platform is not a repository</h2>
            <h3>
                We used GitHub advanced search to find candidate platforms, but
                the unit of analysis is the platform, not the repository
            </h3>
            <div class="slide-body">
                <div class="kpis" style="--n: 3">
                    <div class="panel center">
                        <div class="stat">31</div>
                        <div class="stat-label">
                            repositories in the CODA organisation; none of them
                            a platform-level artefact
                        </div>
                    </div>
                    <div class="panel center">
                        <div class="stat">9</div>
                        <div class="stat-label">
                            repositories for PHT-meDIC / PADME, split across two
                            institutions
                        </div>
                    </div>
                    <div class="panel center">
                        <div class="stat">1</div>
                        <div class="stat-label">
                            for FLIP, Kaapana, vantage6; <br />
                            a monorepo each
                        </div>
                    </div>
                </div>
                <p class="small" style="margin-top: 0.7em">
                    CODA ships as <code>auth-service</code>,
                    <code>fhir-store</code>, <code>hub-api</code>,
                    <code>site-deployer</code> and 27 more; no single one of
                    which is findable as "a federated learning platform".
                </p>
                <p>
                    <strong>Rule:</strong> the unit is the <em>platform</em>.
                    Where a platform spans an organisation, the organisation is
                    the unit and its repositories are read together.
                </p>
            </div>
            <aside class="notes">
                (~1 min 15) This looks like pedantry until you see CODA. It's
                the closest analogue to FLIP in the whole comparison: nine
                Québec hospitals, a million patients, and it is structurally
                invisible to a GitHub repository search. If you don't state the
                unit of analysis, you silently exclude exactly the mature,
                well-engineered platforms that decomposed properly.
            </aside>
        </section>

        <!-- 5 · The shared query design. Stated once here so the four sources are read as one method. -->
        <section>
            <div class="eyebrow">The search · one design</div>
            <h2>One query, three concept blocks, four sources</h2>
            <div class="slide-body">
                <p class="small" style="margin: 0 0 0.4em">
                    The same conjunction is asked of every source. What changes
                    is only what each source
                    <em>forces</em> us to change.
                </p>
                <div class="cols" style="--n: 3">
                    <div class="panel flip">
                        <h3>A · federated computation</h3>
                        <p class="small">
                            federated learning / analytics / evaluation /
                            benchmarking / training / infrastructure ·
                            <em>swarm learning</em> ·
                            <em>personal health train</em> ·
                            <em>distributed analysis</em>
                        </p>
                    </div>
                    <div class="panel">
                        <h3>B · a system, not a study</h3>
                        <p class="small">
                            platform · framework · infrastructure · toolkit ·
                            software · open source · system · tool
                        </p>
                    </div>
                    <div class="panel">
                        <h3>C · health setting</h3>
                        <p class="small">
                            clinical · hospital · patient · medical · health ·
                            biomedical · imaging · radiology · multicentre ·
                            multi-institutional
                        </p>
                    </div>
                </div>
                <p class="center small" style="margin: 0.45em 0 0.35em">
                    <code
                        >A <span class="accent">AND</span> B
                        <span class="accent">AND</span> C</code
                    >, searched over title and abstract
                </p>
                <table style="font-size: 0.44em">
                    <thead>
                        <tr>
                            <th>Source</th>
                            <th>Forced adaptation</th>
                            <th>Why</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong>PubMed</strong></td>
                            <td>
                                reference form; adds the
                                <code>[MeSH]</code> descriptor to A and
                                <code>*</code> wildcards in B and C
                            </td>
                            <td>
                                the only source with a controlled vocabulary
                            </td>
                        </tr>
                        <tr>
                            <td><strong>arXiv</strong></td>
                            <td>
                                A free text only; B and C list singular
                                <em>and</em> plural explicitly
                            </td>
                            <td>no MeSH, and no stemming or wildcards</td>
                        </tr>
                        <tr>
                            <td><strong>medRxiv</strong></td>
                            <td>
                                <strong>C dropped</strong>; A widened to
                                <code>federat*</code> + split learning; matched
                                locally over 109,695 records
                            </td>
                            <td>
                                no search endpoint, and the server is medical by
                                construction
                            </td>
                        </tr>
                        <tr>
                            <td><strong>GitHub</strong></td>
                            <td>
                                not expressible; becomes
                                <strong>18 narrow queries</strong>, unioned, of
                                which 8 name a platform outright
                            </td>
                            <td>
                                five-operator cap; no abstract, only topics and
                                description
                            </td>
                        </tr>
                    </tbody>
                </table>
            </div>
            <aside class="notes">
                (~1 min 30) This slide exists because the four searches should
                look like one method, not four improvisations. Until now the
                deck showed them one at a time, which made them look ad hoc. The
                design is a single A AND B AND C conjunction: is it federated,
                is it a system rather than a study, is it health. Everything
                after this is a deviation forced by a source, and I want each
                one attributable. PubMed is the reference form. arXiv loses MeSH
                and stemming, so the blocks get longer but mean the same thing.
                medRxiv is the interesting one: no search endpoint, so we fetch
                all 109,695 records and match locally, and we drop block C
                because every record on a medical preprint server satisfies it.
                Dropping a block would normally be a red flag, so say plainly
                that it was widened elsewhere to compensate. GitHub can't
                express the conjunction at all. The honest bit is the last row:
                eight queries name a platform outright, which only works for
                platforms you can already name. If someone asks why the
                adaptations differ, the answer is that every one of them is a
                property of the source, not a choice about what we wanted to
                find.
            </aside>
        </section>

        <!-- 6 · The exact queries, as sent. -->
        <section>
            <div class="eyebrow">The search · as sent</div>
            <h2>The four queries, verbatim</h2>
            <div class="slide-body">
                <p class="small" style="margin: 0 0 0.35em">
                    What was actually submitted, so the counts can be
                    reproduced, and argued with. Blocks
                    <span class="accent">A</span> /
                    <span class="accent">B</span> /
                    <span class="accent">C</span> are as defined opposite.
                </p>
                <div class="cols tight" style="--n: 2">
                    <div class="panel flip">
                        <h3>PubMed — 1,490 records</h3>
                        <pre
                            style="
                                font-size: 0.34em;
                                line-height: 1.5;
                                margin: 0.2em 0 0;
                                white-space: pre-wrap;
                            "
                        ><code>("Federated Learning"[MeSH]
  OR "federated learning"[tiab]
  OR "federated analytics"[tiab]
  OR "federated evaluation"[tiab]
  OR "federated benchmarking"[tiab]
  OR "federated training"[tiab]
  OR "federated infrastructure"[tiab]
  OR "swarm learning"[tiab]
  OR "personal health train"[tiab]
  OR "distributed analysis"[tiab])

AND (platform*[tiab] OR framework*[tiab]
  OR infrastructur*[tiab] OR toolkit*[tiab]
  OR software[tiab] OR "open source"[tiab]
  OR open-source[tiab] OR system*[tiab]
  OR tool*[tiab])

AND (clinic*[tiab] OR hospital*[tiab]
  OR patient*[tiab] OR medical[tiab]
  OR health*[tiab] OR biomedical[tiab]
  OR imaging[tiab] OR radiolog*[tiab]
  OR multicent*[tiab] OR multi-cent*[tiab])</code></pre>
                    </div>
                    <div class="panel">
                        <h3>arXiv — 1,190 records</h3>
                        <pre
                            style="
                                font-size: 0.34em;
                                line-height: 1.5;
                                margin: 0.2em 0 0;
                                white-space: pre-wrap;
                            "
                        ><code>((ti:"federated learning" OR abs:"federated learning")
  OR ... 9 phrases, each ti: OR abs:)
AND ((ti:"platform" OR abs:"platform")
  OR ... 14 terms, singular and plural
  spelled out explicitly)
AND ((ti:"clinic" OR abs:"clinic")
  OR ... 15 health terms)

No MeSH, no stemming, no wildcards.
Every surface form is listed.</code></pre>
                        <h3 style="margin-top: 0.5em">
                            medRxiv — 167 of 109,695
                        </h3>
                        <pre
                            style="
                                font-size: 0.34em;
                                line-height: 1.5;
                                margin: 0.2em 0 0;
                                white-space: pre-wrap;
                            "
                        ><code>A  /federat\w*|swarm[- ]learning|
     personal health train|
     distributed analysis|split learning/i

B  /(platforms?|frameworks?|
     infrastructures?|toolkits?|
     software|open[- ]?sourc\w*|
     systems?|tools?)/i

Applied locally: no search endpoint.</code></pre>
                    </div>
                </div>
            </div>
            <aside class="notes">
                (~1 min 30) This is the reproducibility slide, and it is the one
                a methods reviewer will photograph. Everything here was actually
                submitted and the counts are what came back on the run date. Two
                things to say out loud. First, PubMed is shown complete because
                it is short enough to show complete: the [MeSH] descriptor, the
                ten free-text phrases, and the wildcards on blocks B and C.
                Second, arXiv and medRxiv are shown as structure rather than in
                full, because the arXiv query is 1,670 characters of ti:/abs:
                pairs; the ellipses mark where the nine phrases and fourteen
                system terms repeat, and the full strings are in the appendix
                and in search_protocol.sh. Be upfront about that: it is a slide,
                not the supplementary material, and I would rather show an
                honest excerpt with the rule stated than shrink a wall of text
                until nobody can read it. If someone wants the literal arXiv
                string, point them at the appendix; it is reproduced verbatim
                there. The medRxiv row is the interesting contrast: there is no
                query to show, because there is no search endpoint, so what we
                publish is the matching rule and the fact that it was applied to
                all 109,695 records.
            </aside>
        </section>

        <!-- 7 · Inclusion rules. -->
        <section>
            <div class="eyebrow">Criteria · 2 of 4</div>
            <h2>Inclusion — all five must hold</h2>
            <div class="slide-body">
                <ol class="contribs tight">
                    <li>
                        A <strong>named, reusable software system</strong>:
                        platform, infrastructure or FL engine. Not a study, an
                        algorithm or a dataset.
                    </li>
                    <li>
                        Computation distributed to
                        <strong>two or more independently governed sites</strong
                        >, with record-level data remaining at each.
                    </li>
                    <li>
                        <strong>Health-designed or health-evaluated</strong>:
                        <em>or</em> a general-purpose engine documented as
                        consumed by a health platform in the set.
                    </li>
                    <li>
                        Described by a
                        <strong
                            >peer-reviewed publication, a public repository, or
                            both</strong
                        >.
                    </li>
                    <li>
                        Publicly available at some point between
                        <strong>2014 and 2026</strong>.
                    </li>
                </ol>
                <div class="flip-note">
                    <span class="flip-tag">Note</span> Rule 3's second clause
                    admits the engine layer: Flower, FLARE, OpenFL. FLIP
                    <em>consumes</em> those rather than competing with them.
                </div>
            </div>
            <aside class="notes">
                (~1 min 30) Walk these. The one to dwell on is rule 3; without
                its second clause the engines fall out, and then the table can't
                show that we're a consumer at that layer rather than a rival.
                2014 as the start year isn't arbitrary: it's DataSHIELD, which
                set the disclosure control discipline this whole field
                inherited.
            </aside>
        </section>

        <!-- 8 · Exclusion. -->
        <section>
            <div class="eyebrow">Criteria · 3 of 4</div>
            <h2>Exclusion</h2>
            <div class="slide-body">
                <ul class="crosslist">
                    <li>
                        Method papers with <strong>no released system</strong>.
                    </li>
                    <li>
                        <strong>Simulation-only frameworks</strong> with no site
                        component anyone could install.
                    </li>
                    <li>Surveys, reviews, position papers, protocols.</li>
                    <li>
                        <strong>Single-institution portals</strong> with no
                        federated computation.
                    </li>
                    <li>
                        Systems
                        <strong
                            >superseded by a successor already in the
                            set</strong
                        >
                        we score the successor.
                    </li>
                </ul>
                <p class="small muted" style="margin-top: 0.6em">
                    The last rule is why the table carries
                    <strong>NEBULA</strong> rather than Fedstellar, one project
                    under a new name rather than two competitors, and why it
                    carries a single <strong>PHT-meDIC / PADME</strong> row: two
                    distinct codebases sharing a star topology and a
                    station-approval model.
                </p>
            </div>
            <aside class="notes">
                (~1 min) Faster through these. The supersession rule earns its
                place: Fedstellar was renamed NEBULA, and if you don't state the
                rule you look like you cited a paper whose software no longer
                exists under that name. Carry both names so the citation still
                resolves, and so the reader can find the work under either
                title. Note the two deliberate exclusions that are not
                judgements of quality: FEM-orchestrator moved its source to
                GitLab and monai-fl is an aspirational example, so neither says
                anything about whether the software is good, only about whether
                we could score it.
            </aside>
        </section>

        <!-- 9 · Evidence tiers — this is what makes the star marker principled. -->
        <section>
            <div class="eyebrow">Criteria · 4 of 4</div>
            <h2>
                Evidence tiers, and what
                <span class="accent">⋆</span> means
            </h2>
            <div class="slide-body">
                <div class="cols" style="--n: 3">
                    <div class="panel flip">
                        <h3>Tier A · scored from code</h3>
                        <p class="small">
                            Source we could read. Every mark traces to a file we
                            can name.
                        </p>
                    </div>
                    <div class="panel fla3">
                        <h3>
                            Tier B · scored from paper
                            <span class="accent">⋆</span>
                        </h3>
                        <p class="small">
                            Described in a publication, but code closed,
                            unreleased, or hosted where we could not reach it.
                            Marked, and read as the <em>authors' claim</em>.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Tier C · excluded</h3>
                        <p class="small">
                            Neither a publication nor reachable code.
                        </p>
                    </div>
                </div>
                <p class="small" style="margin-top: 0.7em">
                    Inclusion rule 4 is deliberately <em>disjunctive</em>.
                    Scoring Tier B from code is impossible; excluding it would
                    silently drop the commercial and consortium platforms a
                    hospital is most likely to be offered.
                </p>
                <div class="flip-note">
                    <span class="flip-tag">Result</span> Scoring from local
                    clones leaves Tier B at <strong>three rows</strong>:
                    GenoMed4All, INCISIVE and TrainTracks. Everything else is
                    read from source at a recorded commit.
                </div>
            </div>
            <aside class="notes">
                (~1 min 15) This slide retro-fits a principle onto a marker we
                were already using. Before this, the star was "we couldn't
                check"; an apology. Now it's a tier in a stated scheme, which is
                a much stronger position: we're not hiding the weakness, we're
                classifying it. Two honest edges to state. FeatureCloud's paper
                describes its architecture as open source and we could not
                locate frontend, backend or controller, so its web UI cell is a
                dash rather than a score. Say that as a limit on what we could
                verify, not as a claim about their code. And Tier B rows are
                read as the authors' claim, which is a weaker form of evidence,
                but excluding them would hide the commercial and consortium
                options a hospital is most likely to be offered. The box is the
                good news: only three rows in the whole table rest on the
                authors' claim rather than on source we read.
            </aside>
        </section>

        <!-- 10 · PubMed query. -->
        <section>
            <div class="eyebrow">Search · source 1</div>
            <h2>PubMed — the reference form</h2>
            <div class="slide-body">
                <p class="small">
                    The three blocks as submitted, with the one addition no
                    other source has: the <code>[MeSH]</code> descriptor in
                    block A.
                </p>
                <div class="kpis" style="--n: 1; margin-top: 0.5em">
                    <div class="panel center">
                        <div class="stat">1,485</div>
                        <div class="stat-label">
                            records returned by A
                            <span class="accent">AND</span> B
                            <span class="accent">AND</span> C
                        </div>
                    </div>
                </div>
                <p class="small muted" style="margin-top: 0.6em">
                    Recovers <strong>13 in the field</strong>: APPFL (as
                    APPFLx), CODA, FeTS, INCISIVE, PHT&minus;meDIC/PADME,
                    MedPerf, FLARE, OpenFL, FedKBP<sup>+</sup>, vantage6,
                    DataSHIELD, gLinDA and TrainTracks.
                </p>
            </div>
            <aside class="notes">
                (~1 min) We've already seen the block design, so this slide is
                just the count and what it caught. The interesting part is what
                it took to get block A right. Next slide. Note we release the
                executable query, so anyone can re-run this and get a number to
                argue with.
            </aside>
        </section>

        <!-- 12 · GitHub, and the operator ceiling. -->
        <section>
            <div class="eyebrow">Search · source 2</div>
            <h2>GitHub — a union, because one query is not allowed</h2>
            <div class="slide-body">
                <p class="small">
                    The search API rejects any query using
                    <strong>more than five boolean operators</strong>. So the
                    strategy is ten narrow queries, unioned, not one expression.
                </p>
                <table style="font-size: 0.5em">
                    <thead>
                        <tr>
                            <th>Query</th>
                            <th>Catches</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>
                                <code
                                    >topic:federated-learning …
                                    pushed:&gt;=2024-09-01 stars:&gt;=20</code
                                >
                            </td>
                            <td>the mainstream of the field</td>
                        </tr>
                        <tr>
                            <td>
                                <code>topic:federated-learning-framework</code>
                            </td>
                            <td>APPFL, which does not carry the plain topic</td>
                        </tr>
                        <tr>
                            <td>
                                <code>topic:federated topic:healthcare</code>
                            </td>
                            <td>
                                Fed-BioMed, whose topics are split into single
                                words
                            </td>
                        </tr>
                        <tr>
                            <td>
                                <code
                                    >topic:federated-learning
                                    topic:medical-imaging</code
                                >
                            </td>
                            <td>the imaging platforms</td>
                        </tr>
                        <tr>
                            <td>
                                <code
                                    >"federated learning"
                                    in:name,description</code
                                >
                            </td>
                            <td>
                                Kaapana and MedPerf, which carry no useful
                                topics but a real description
                            </td>
                        </tr>
                        <tr>
                            <td>
                                <code>"federated evaluation"</code> ·
                                <code>"personal health train"</code> ·
                                <code>"swarm learning"</code>
                            </td>
                            <td>
                                the sub-literatures with their own vocabulary
                            </td>
                        </tr>
                        <tr>
                            <td>
                                <code>repo:vantage6/vantage6</code> ·
                                <code>repo:datashield/dsBase</code> · six more
                            </td>
                            <td>
                                the metadata-poor ones, reachable
                                <strong>only</strong> by name
                            </td>
                        </tr>
                    </tbody>
                </table>
                <div class="kpis" style="--n: 2; margin-top: 0.5em">
                    <div class="panel center">
                        <div class="stat">18</div>
                        <div class="stat-label">queries</div>
                    </div>
                    <div class="panel center">
                        <div class="stat">2,330</div>
                        <div class="stat-label">
                            unique repositories in the union
                        </div>
                    </div>
                </div>
                <p class="small muted" style="margin-top: 0.6em">
                    Recovers <strong>20 of the 22 with a public repo</strong>.
                    The last eight queries target platforms <em>by name</em>:
                    which is what recovered FLA<sup>3</sup>, FeatureCloud,
                    vantage6 and DataSHIELD, none of which carry a matching
                    topic or description. The two still missed are
                    FedKBP<sup>+</sup>, an empty placeholder repository, and
                    gLinDA, hosted on a university GitLab.
                </p>
            </div>
            <aside class="notes">
                (~1 min 15) The five-operator cap is a real constraint people
                don't know about until they hit a 422. Point at the middle rows
                because every one of them exists because a specific platform was
                being missed, and the reason was always maintainer metadata
                rather than anything about the software. Fed-BioMed's topics are
                literally "federated" and "learning" as two separate tags. The
                honest limit to state: topic and description queries alone
                cannot reach a repository that fills in neither, at any
                threshold. That is why eight of the queries target names, and it
                is the difference between recall of 16 and 20.
            </aside>
        </section>

        <!-- 13 · Sources 3 & 4: the preprint layer. -->
        <section>
            <div class="eyebrow">Search · sources 3 &amp; 4</div>
            <h2>arXiv and medRxiv — the preprint layer</h2>
            <div class="slide-body">
                <div class="cols" style="--n: 2">
                    <div class="panel">
                        <h3>arXiv · 1,190 records</h3>
                        <p class="small">
                            Where the FL <strong>engines</strong> publish before
                            a peer-reviewed venue. A preprint counts as evidence
                            only when it also names a repo or a journal.
                        </p>
                        <p class="small muted">
                            <strong>Recovers in the field:</strong>
                            FeatureCloud, FLA<sup>3</sup>,
                            PHT&minus;meDIC/PADME, MedPerf, FeTS, DataSHIELD,
                            Swarm Learning, FLARE, OpenFL, Flower, Fed-BioMed,
                            Substra, FedML and FedKBP<sup>+</sup>: 14 of 30.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>medRxiv · 167 matches</h3>
                        <p class="small">
                            All <strong>109,695</strong> records since the 2019
                            launch, matched <strong>locally</strong>, because
                            there is no search endpoint to submit a query to.
                        </p>
                        <p class="small muted">
                            <strong>Recovers in the field:</strong> GenoMed4All
                            invisible to the other three because its only
                            description is a preprint, plus
                            PHT&minus;meDIC/PADME, DataSHIELD, Swarm Learning
                            and vantage6: 5 of 30.
                        </p>
                    </div>
                </div>
                <div class="kpis" style="--n: 4; margin-top: 0.8em">
                    <div class="panel center">
                        <div class="stat">13</div>
                        <div class="stat-label">PubMed / 30</div>
                    </div>
                    <div class="panel center">
                        <div class="stat">20</div>
                        <div class="stat-label">GitHub / 22 w. repo</div>
                    </div>
                    <div class="panel center">
                        <div class="stat">14</div>
                        <div class="stat-label">arXiv / 30</div>
                    </div>
                    <div class="panel center">
                        <div class="stat">5</div>
                        <div class="stat-label">medRxiv / 30</div>
                    </div>
                </div>
                <p style="margin-top: 0.5em">
                    Each source sees a different slice of the field:
                    <strong>together the four recover all 30</strong>.
                </p>
            </div>
            <aside class="notes">
                (~45 s) The engines are on arXiv, the very-new clinical
                platforms are on medRxiv, and neither index accepts the other's
                content. This is the slide that justifies four sources instead
                of two: each one pays for itself in exactly one platform the
                others cannot see.
            </aside>
        </section>

        <!-- 14 · The scoreboard. The core empirical result of the methods work. -->
        <section>
            <div class="eyebrow">Validation</div>
            <h2>Does either search actually find them?</h2>
            <div class="slide-body">
                <table style="font-size: 0.5em">
                    <thead>
                        <tr>
                            <th>Source</th>
                            <th>Pool</th>
                            <th>Recall</th>
                            <th>Blind spot</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong>PubMed</strong></td>
                            <td>1,485 records</td>
                            <td><span class="cov partial">◐</span>13 / 30</td>
                            <td>FL engines: <strong>2 of 8</strong></td>
                        </tr>
                        <tr>
                            <td><strong>GitHub</strong></td>
                            <td>2,330 repositories</td>
                            <td>
                                <span class="cov partial">◐</span>20 / 22 with
                                repos
                            </td>
                            <td>empty placeholders; non-GitHub hosting</td>
                        </tr>
                        <tr>
                            <td><strong>arXiv</strong></td>
                            <td>1,190 preprints</td>
                            <td><span class="cov partial">◐</span>14 / 30</td>
                            <td>health-restricted; misses engines</td>
                        </tr>
                        <tr>
                            <td><strong>medRxiv</strong></td>
                            <td>167 matches</td>
                            <td><span class="cov partial">◐</span>5 / 30</td>
                            <td>
                                catches GenoMed4All, invisible to the other
                                three
                            </td>
                        </tr>
                        <tr class="is-ours">
                            <td><strong>All four</strong></td>
                            <td>—</td>
                            <td>
                                <span class="cov yes">●</span
                                ><strong>30 / 30</strong>
                            </td>
                            <td>nothing missed</td>
                        </tr>
                    </tbody>
                </table>
                <p style="margin-top: 0.7em">
                    No single source recovers even two-thirds of the field.
                    <strong>They fail in different directions</strong>, and it
                    takes all four to reconstruct the full set. The last gap,
                    GenoMed4All, is a medRxiv preprint that PubMed, GitHub and
                    arXiv cannot see.
                </p>
            </div>
            <aside class="notes">
                (~1 min 30) The headline. Pause on the middle column. If we had
                run a PubMed-only search and called it systematic, we'd have
                published a comparison missing five of the eight FL engines,
                including Substra and Fed-BioMed. GitHub's blind spot is
                structural, not a matter of tuning: metadata-poor clinical
                platforms are reached by the name-targeted queries, and what
                remains beyond reach is a repository with no code in it and a
                platform that isn't on GitHub at all. The last row is the good
                news, but it is a union of four sources, and I'll be explicit
                about that in three slides.
            </aside>
        </section>

        <!-- 16 · The popularity trap, quantified. -->
        <section>
            <div class="eyebrow">Validation</div>
            <h2>The popularity trap</h2>
            <div class="slide-body">
                <p class="small">
                    Median GitHub stars, by layer of the comparison:
                </p>
                <div class="kpis" style="--n: 3">
                    <div class="panel center">
                        <div class="stat">903</div>
                        <div class="stat-label">FL engines</div>
                    </div>
                    <div class="panel center">
                        <div class="stat cyan">110</div>
                        <div class="stat-label">
                            deployed healthcare platforms
                        </div>
                    </div>
                    <div class="panel center">
                        <div class="stat pink">32</div>
                        <div class="stat-label">
                            federated analytics infrastructures
                        </div>
                    </div>
                </div>
                <p style="margin-top: 0.8em">
                    Any star floor tight enough to make the candidate pool
                    screenable
                    <strong
                        >removes the clinically deployed platforms first</strong
                    >
                    exactly the population the paper is about. A
                    <code>stars:&gt;=100</code> floor recovers
                    <strong>6 of 18</strong>.
                </p>
                <div class="flip-note">
                    <span class="flip-tag">But</span> Dropping the floor
                    entirely grew the pool 579 &rarr; 2,107 and recovered
                    <strong>nothing</strong>. The star floor was never the
                    binding constraint; the queries only ever read topics and
                    two generic name phrases. The floor is now
                    <strong>&ge;1 contributor</strong>, applied after the fetch
                    as a liveness check: GitHub has no collaborator qualifier,
                    so it cannot live in the query.
                </div>
            </div>
            <aside class="notes">
                (~1 min 30) This is the one I'd generalise beyond our paper.
                Popularity filtering is the default hygiene step in any
                repository-mining study, and in this domain it is actively
                biased: hospital software is written by small teams for
                institutional users who don't star repositories. The 6-of-18
                number is measured, not hypothetical. But the twist in the box
                is the more useful lesson, and it cost me a day: I assumed the
                floor was the problem, removed it, quadrupled the pool and
                recovered not one missing platform. The instrument was blind in
                a different dimension than the one I was tuning. Measure which
                constraint actually binds before you relax it. If anyone asks
                why the floor is &ge;1 and not &ge;2: two contributors would
                drop FLA<sup>3</sup>, which is single-author.
            </aside>
        </section>

        <!-- 17 · The honest conclusion about the protocol. -->
        <section>
            <div class="eyebrow">Validation</div>
            <h2>Hybrid protocol</h2>
            <div class="slide-body">
                <ol class="contribs tight">
                    <li>
                        <strong>Query-driven identification</strong> from both
                        sources: reproducible, bounded, released as a script.
                    </li>
                    <li>
                        <strong>Citation tracking</strong> over included records
                        and recent surveys. It completes the metadata-poor and
                        CS-venue remainder.
                    </li>
                    <li>
                        <strong>Eligibility criteria</strong> applied to the
                        union, with evidence tiers.
                    </li>
                    <li>
                        <strong>Report the recall of each channel</strong>, so
                        the reader knows what the instrument can and cannot see.
                    </li>
                </ol>
            </div>
            <aside class="notes">
                (~1 min) The framing that matters: all four pools are screened.
                PubMed 1,483/1,485; arXiv 1,189/1,190; medRxiv 126/126&mdash;
                each record was judged against the five inclusion criteria and a
                per-record verdict file written. If anyone spots that the recall
                table says 167 medRxiv matches and this says 126: the screened
                pool is the 5 September run, and a re-run on 10 September
                returned 167. The difference is new deposits, not a change of
                query. GitHub is screened on the search metadata (title,
                description, language, topics) for all 2,330 in the pool, 1,773
                of those from README plus metadata; 31 repos were cloned and
                read from source. What we are not claiming is a
                systematic-review-grade pass over the literature&mdash; that is
                a validated, reproducible identification protocol with a
                measured recall per channel, and a citation-tracking pass on
                top. That is a stronger and more honest position than a
                screening claim we cannot stand behind.
            </aside>
        </section>

        <!-- 18 · The nine columns, before the tables that use them. -->
        <section>
            <div class="eyebrow">The comparison</div>
            <h2>Nine columns, each with a falsifiable test</h2>
            <div class="slide-body">
                <table style="font-size: 0.5em">
                    <tbody>
                        <tr>
                            <td><strong>Site data model</strong></td>
                            <td>
                                the schema a site must present; imposing none is
                                <em>Flexible</em>
                            </td>
                        </tr>
                        <tr>
                            <td><strong>Peer-to-peer</strong></td>
                            <td>
                                full decentralisation, no central coordinator
                            </td>
                        </tr>
                        <tr>
                            <td><strong>Web UI runs the workflow</strong></td>
                            <td>
                                a browser executes a <em>complete</em> project:
                                not monitoring or admin only
                            </td>
                        </tr>
                        <tr>
                            <td><strong>Pluggable FL engine</strong></td>
                            <td>
                                &gt;1 FL <em>engine</em> behind one contract.
                                Supporting many ML frameworks does not count
                            </td>
                        </tr>
                        <tr>
                            <td><strong>Federated cohort query</strong></td>
                            <td>
                                query executed at each site, aggregate-only
                                return
                            </td>
                        </tr>
                        <tr>
                            <td><strong>On-demand PACS imaging</strong></td>
                            <td>
                                DICOM retrieval from an institutional PACS; bulk
                                pre-ingestion scores partial
                            </td>
                        </tr>
                        <tr>
                            <td><strong>Per-site project approval</strong></td>
                            <td>
                                per-project, per-site approval recorded as an
                                object that <em>gates execution</em>
                            </td>
                        </tr>
                        <tr>
                            <td><strong>Runtime policy enforced</strong></td>
                            <td>code that blocks, not a human sign-off</td>
                        </tr>
                        <tr>
                            <td><strong>Exportable provenance</strong></td>
                            <td>
                                structured end-to-end record for a completed
                                job; logs and dashboards don't qualify
                            </td>
                        </tr>
                    </tbody>
                </table>
                <p class="small muted" style="margin-top: 0.5em">
                    Every one is phrased so a reader with the source could prove
                    us wrong.
                </p>
            </div>
            <aside class="notes">
                (~1 min 30) Don't read all nine; pick three. "Pluggable FL
                engine" is worth a beat: a co-author scored FLARE, Flower and
                OpenFL as full marks there, reading it as "supports many ML
                frameworks". They were right under that reading, which told us
                the column name was the defect. We renamed it and wrote the
                exclusion into the criteria. If a co-author misreads a header,
                reviewers will too.
            </aside>
        </section>

        <!-- 19 · Capability table, part 1: our own peer group. Ten rows — the
             systems a hospital would actually weigh against FLIP. Five rows
             that used to sit here are now on slide 20. Four are in the
             "General-purpose FL frameworks" band, because each one's own
             manifest describes a library or framework rather than a deployment:
             APPFL is "an open-source and highly extensible software framework",
             FL4Health "contains the federated learning (FL) engine" and is a
             library, FeatureCloud ships a pip app-engine with no self-hostable
             frontend, backend or controller, and Starfish-FL is "an agentic
             federated learning framework" aimed at healthcare, resource
             allocation and finance alike. None of the four exposes a clinical
             connector (no PACS, FHIR or OMOP surface), which is what kept them
             out of a hospital's shortlist. Keeping them in the set matters: the
             searches recovered them, and rule 3 admits them as health-designed
             or health-evaluated. The fifth is MEDIATA, now in the "Federated
             analytics infrastructures" band: it harmonises and queries data
             across sites, and there is no model artefact or model registry
             anywhere in its three repositories, so it is a neighbour on the
             cohort-query column rather than a competitor at the training layer.

             Licence column (both capability tables now carry one): read from
             the repository itself at the commit the platform was scored
             against, never from memory or a project's marketing page. The red
             marker is `cov no`, the same one the capability columns use,
             meaning a hospital cannot adopt the row as-is — copyleft, a
             non-OSI research licence, or no licence at all. Where no licence
             exists the cell reads Copyright, because the Berne default is all
             rights reserved: "we found no source" and "there is no permission
             to use this" are the same fact seen from a hospital's side. Rows
             needing a note:
               · CODA is 32 repositories under one org; 21 carry GPL-3.0 and
                 none carries anything else, so the row reads GPL-3.0.
               · FeTS ships a composite licence — a BSD-style UI licence over an
                 Apache-2.0 federated backend — hence "BSD + Apache".
               · GenoMed4All and INCISIVE publish neither a repository nor a
                 licence, so their code is all rights reserved by default; the
                 cell reads Copyright rather than a dash, because "we found no
                 source" and "there is no permission to use this" are the same
                 fact seen from a hospital's side.
             FL4Health's Vector Institute licence (academic, sponsor or partner
             use only — not OSI-approved) travelled with it to slide 20.
             Provenance for slide 20's column, which was read the same way:
             Swarm Learning keeps its Apache-2.0 text under docs/Generic/
             LICENSE.md rather than at the repository root, and additionally
             requires HPE's AutoPass licence server to run (a free community
             edition exists, but it needs an HPE account); gLinDA lives on a
             university GitLab rather than GitHub and declares BSD-3 in its
             setup.py; FEDn's Apache-2.0 covers the client SDK only, as the
             combiner and controller tiers now ship as commercial images;
             FedKBP+ and TheODen carry no licence file and GitHub detects none,
             so they read Copyright. -->

        <section>
            <div class="eyebrow">
                The comparison · deployed healthcare platforms
            </div>
            <h2>Our peer group</h2>
            <div class="slide-body">
                <p class="small" style="margin: 0 0 0.3em">
                    Deployed healthcare platforms: the closest neighbours, and
                    the rows a hospital would actually weigh against FLIP.
                </p>
                <table class="cap-table" style="font-size: 0.44em">
                    <colgroup>
                        <col style="width: 15.5%" />
                        <col style="width: 13.5%" />
                        <col style="width: 11%" />
                        <col span="8" style="width: 7.5%" />
                    </colgroup>
                    <thead>
                        <tr>
                            <th>Platform</th>
                            <th>Data model</th>
                            <th>Licence</th>
                            <th>P2P</th>
                            <th>Web UI</th>
                            <th>Engines</th>
                            <th>Cohort</th>
                            <th>PACS</th>
                            <th>Approval</th>
                            <th>Policy</th>
                            <th>Prov.</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr class="is-ours">
                            <td>
                                <a
                                    href="https://github.com/londonaicentre/FLIP"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    ><strong>FLIP (ours)</strong></a
                                >
                            </td>
                            <td>OMOP + MI-CDM</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/kaapana/kaapana"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >JIP / Kaapana</a
                                >
                            </td>
                            <td>DICOM</td>
                            <td><span class="cov no" title="AGPL-3.0 — strong network copyleft: deploying a modified version as a service obliges you to publish your source">AGPL&#8209;3.0</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/mlcommons/medperf"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >MedPerf</a
                                >
                            </td>
                            <td>Benchmark-defined</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/PHT-Medic/central"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >PHT-meDIC / PADME</a
                                >
                            </td>
                            <td>Flexible (FHIR payload)</td>
                            <td>MIT</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/coda-platform"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >CODA</a
                                >
                            </td>
                            <td>FHIR + DICOM</td>
                            <td><span class="cov no" title="GPL-3.0 — copyleft: derivative works must be released under the same licence">GPL&#8209;3.0</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/bloodcounts/FLAAA"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >FLA<sup>3</sup></a
                                >
                            </td>
                            <td>Tabular</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://doi.org/10.1101/2025.08.07.25333044"
                                    target="_blank"
                                    rel="noopener"
                                    title="Article"
                                    >GenoMed4All<sup class="accent">⋆</sup></a
                                >
                            </td>
                            <td>FHIR</td>
                            <td><span class="cov no" title="No licence located for this project — under the Berne Convention the default is all rights reserved">Copyright</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/NeuroFlame/NeuroFLAME"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >NeuroFLAME</a
                                >
                            </td>
                            <td>NIfTI + CSV</td>
                            <td>MIT</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/FeTS-AI/Front-End"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >FeTS</a
                                >
                            </td>
                            <td>NIfTI / BraTS</td>
                            <td>BSD + Apache</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://cordis.europa.eu/project/id/952179"
                                    target="_blank"
                                    rel="noopener"
                                    title="Article"
                                    >INCISIVE<sup class="accent">⋆</sup></a
                                >
                            </td>
                            <td>bespoke CDM</td>
                            <td><span class="cov no" title="No licence located for this project — under the Berne Convention the default is all rights reserved">Copyright</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                    </tbody>
                </table>
                <p class="small muted" style="margin-top: 0.4em">
                    <span class="cov yes">●</span>ships &nbsp;
                    <span class="cov partial">◐</span>partial or by extension
                    &nbsp; <span class="cov no">✕</span>not supported &nbsp;
                    <span class="muted">—</span> no reachable source settles it
                    &nbsp; <span class="muted">n/a</span> not applicable
                </p>
            </div>
            <aside class="notes">
                (~2 min) The densest slide; give it time, and mind the tone;
                these are other teams' shipped products, several in clinical
                use. Three things to point at. One — the cohort-query and PACS
                columns are nearly empty, and that combination is where FLIP
                fits. CODA is the only other full mark on cohort query, and it
                went to FHIR where we went to OMOP. MEDIATA reaches the same
                column from the analytics side and now sits in that band on the
                next table, with no training engine at all. Both are solving it
                well for their own question. Two — Kaapana scores full on PACS and on
                per-site approval, and its runtime policy is stronger than ours:
                Traefik forwardAuth against an Open Policy Agent decision,
                blocking the request path, not the rendered page. That is a
                design to learn from, not a gap to paper over. Three — scoring
                every cell from source means four interface claims in the
                literature get a hard score instead, and three of those score
                lower than the paper implies (FLA³ and FeTS have no web
                workflow here; APPFL and FeatureCloud are the other two, and
                both are now in the frameworks band on the next table). Say
                that plainly: the rule cuts against our own row
                too, and we applied it the same way everywhere. FLA³ is worth a
                positive note here: a real XACML policy engine is the strongest
                governance score in the table.
            </aside>
        </section>

        <!-- 20 · Capability table, part 2. -->
        <section>
            <div class="eyebrow">The comparison · other layers</div>
            <h2>Analytics, peer-to-peer, frameworks and the engines</h2>
            <div class="slide-body">
                <table class="cap-table" style="font-size: 0.40em">
                    <colgroup>
                        <col style="width: 15.5%" />
                        <col style="width: 13.5%" />
                        <col style="width: 11%" />
                        <col span="8" style="width: 7.5%" />
                    </colgroup>
                    <thead>
                        <tr>
                            <th>Platform</th>
                            <th>Data model</th>
                            <th>Licence</th>
                            <th>P2P</th>
                            <th>Web UI</th>
                            <th>Engines</th>
                            <th>Cohort</th>
                            <th>PACS</th>
                            <th>Approval</th>
                            <th>Policy</th>
                            <th>Prov.</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td colspan="10" class="muted">
                                <em>Decentralised (peer-to-peer)</em>
                            </td>
                        </tr>
                        <tr>
                            <td>Swarm Learning</td>
                            <td>Flexible</td>
                            <td><span title="Code is Apache-2.0. Running the framework additionally requires HPE's AutoPass licence server; a free community edition exists but needs an HPE Passport account">Apache&#8209;2.0</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>NEBULA / Fedstellar</td>
                            <td>Built-in only</td>
                            <td><span class="cov no" title="AGPL-3.0 — strong network copyleft: running a modified version as a service obliges you to publish your source">AGPL&#8209;3.0</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>FedKBP<sup>+</sup></td>
                            <td class="muted">—</td>
                            <td><span class="cov no" title="No licence file in the repository and none detected — all rights reserved by default">Copyright</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td>TheODen</td>
                            <td>Flexible</td>
                            <td><span class="cov no" title="No licence file in the repository and none detected — all rights reserved by default">Copyright</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td colspan="10" class="muted">
                                <em>Federated analytics infrastructures</em>
                            </td>
                        </tr>
                        <tr>
                            <td>vantage6</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>TrainTracks<sup class="accent">⋆</sup></td>
                            <td>DataLad / BIDS</td>
                            <td><span class="cov no" title="No repository and no licence located — the row is publication-only (⋆), so there is no published permission to use the code">Copyright</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>DataSHIELD</td>
                            <td>Opal / Armadillo</td>
                            <td><span class="cov no" title="GPL-3.0 — copyleft: derivative works must be released under the same licence">GPL&#8209;3.0</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">—</td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td>gLinDA</td>
                            <td>count matrix</td>
                            <td>BSD&#8209;3</td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td>
                                <a
                                    href="https://github.com/tecnomod-um/MEDIATA_project"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >MEDIATA</a
                                >
                            </td>
                            <td>Tabular + DCAT</td>
                            <td>MIT</td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">—</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td colspan="10" class="muted">
                                <em
                                    >General-purpose FL frameworks:
                                    health-designed or health-evaluated, not
                                    clinical deployments</em
                                >
                            </td>
                        </tr>
                        <tr>
                            <td><a
                                    href="https://github.com/denoslab/starfish-fl"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >Starfish-FL</a
                                ></td>
                            <td>CSV / image ZIP</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td><a
                                    href="https://github.com/APPFL/APPFL"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >APPFL</a
                                ></td>
                            <td class="cov no">Flexible</td>
                            <td>MIT</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td><a
                                    href="https://github.com/FeatureCloud/FeatureCloud"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >FeatureCloud</a
                                ></td>
                            <td class="cov no">Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">—</td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td><a
                                    href="https://github.com/VectorInstitute/FL4Health"
                                    target="_blank"
                                    rel="noopener"
                                    title="Repository"
                                    >FL4Health</a
                                ></td>
                            <td class="cov no">Flexible</td>
                            <td><span class="cov no" title="Vector Institute licence — academic, sponsor or partner use only; not an OSI-approved open-source licence">Vector Inst.</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                        <tr>
                            <td colspan="10" class="muted">
                                <em
                                    >FL engines: FLIP consumes rather than
                                    competes with these</em
                                >
                            </td>
                        </tr>
                        <tr>
                            <td>FLARE</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>Flower</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                        </tr>
                        <tr>
                            <td>OpenFL</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>Fed-BioMed</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>Substra</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">—</td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                        </tr>
                        <tr>
                            <td>FEDn</td>
                            <td>Flexible</td>
                            <td><span title="The client SDK is Apache-2.0; the combiner and controller server tier ships as commercial images, so the full stack cannot be self-hosted from source">Apache&#8209;2.0</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov partial">◐</span></td>
                        </tr>
                        <tr>
                            <td>FATE</td>
                            <td>tabular</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov yes">●</span></td>
                            <td><span class="cov yes">●</span></td>
                        </tr>
                        <tr>
                            <td>FedML</td>
                            <td>Flexible</td>
                            <td>Apache&#8209;2.0</td>
                            <td><span class="cov partial">◐</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td class="muted center">n/a</td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                            <td><span class="cov no">✕</span></td>
                        </tr>
                    </tbody>
                </table>
            </div>
            <aside class="notes">
                (~2 min) Faster. Worth a beat each. Take the new band first:
                four rows moved out of the peer group because they are
                frameworks, not deployments. APPFL and FL4Health are libraries
                — FL4Health's own README calls its contents "the federated
                learning (FL) engine" — FeatureCloud is a pip app-engine whose
                controller we could not find source for, and Starfish-FL is the
                STARFISH project's framework for healthcare, resource allocation
                and finance at once. None of the four speaks to a PACS or a FHIR
                store, which is exactly why a hospital never weighs them against
                FLIP. Two blanks to name while you are here: FeatureCloud's web
                UI stays a dash because no source settles it, and FL4Health's
                licence is academic-only, so it is the one row on either table a
                hospital could not adopt at all. One row joined the analytics
                band for the same reason: MEDIATA, which harmonises and queries
                data across sites but has no model artefact or registry anywhere
                in its repositories — it reaches the cohort-query column without
                ever training a model, which is the analytics layer's signature
                rather than a training platform's. Then the licence column, new
                on this table too. Red means a hospital cannot adopt the row
                as-is, and the ones here are NEBULA's AGPL-3.0, DataSHIELD's
                GPL-3.0, FL4Health's academic-only licence, and Copyright on the
                three rows publishing no licence at all: TrainTracks, which has
                no repository, plus FedKBP+ and TheODen. Worth saying out loud
                that being out of commercial reach says nothing about the
                quality of the research. Note also that Swarm Learning's
                Apache-2.0 sits under docs/Generic/ rather than at the root, and
                that running it needs HPE's licence server. The engines column
                is shaded n/a for engines; you can't ask an engine whether it
                hosts multiple engines. FLARE actually does, via app_opt/flower, which
                wraps a Flower application in FLARE's own FedJob; the shading
                reflects how we partitioned the field, not a limit in FLARE.
                Credit where it is due on governance: Flower ships a typed
                LogEntry audit model with actor, event and status; FLARE's
                ComponentPathAuthorizer refuses to build an unsafe component
                config; Fed-BioMed and FATE both gate a run node-side before any
                work starts; DataSHIELD enforces disclosure control in the
                analysis path, which is the discipline this whole field
                inherited. Substra's policy score drops to partial now that the
                ledger is gone, but its provenance model is still among the most
                complete in the set. The point of this slide is not that engines
                are weak. It is that they stop at governance by design, because
                that is the layer a platform like FLIP sits on.
            </aside>
        </section>

        <!-- 21 · The scoring rule, and what it caught. -->
        <section>
            <div class="eyebrow">Scoring</div>
            <h2>Scored on code, not on the paper</h2>
            <div class="slide-body">
                <p class="small">
                    A capability that exists in a publication but not in a
                    repository cannot be deployed at a hospital. Every cell is
                    scored from a local clone at a recorded commit:
                    <strong>200 cells</strong> across the table. Four where the
                    published description and the source flatly disagree:
                </p>
                <table style="font-size: 0.5em">
                    <thead>
                        <tr>
                            <th>Platform</th>
                            <th>Published as</th>
                            <th>What the source shows today</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong>Substra</strong></td>
                            <td>
                                traceable framework, operations on a distributed
                                ledger
                            </td>
                            <td>
                                orchestrator 0.38.0:
                                <em
                                    >"BREAKING: remove all code related to
                                    <code>distributed</code> mode"</em
                                >. Lineage survives in one central Postgres,
                                with no longer tamper-evident
                            </td>
                        </tr>
                        <tr>
                            <td><strong>FEDn</strong></td>
                            <td>
                                combiner and controller tiers a site operates
                                itself
                            </td>
                            <td>
                                reset to a client SDK, renamed
                                <code>scaleout-client</code>; server tier now
                                commercial images. We score release 0.33.0 and
                                say so
                            </td>
                        </tr>
                        <tr>
                            <td><strong>OpenFL</strong></td>
                            <td>
                                Linux Foundation library for federated training
                                and validation
                            </td>
                            <td>
                                README:
                                <em
                                    >"no longer under active development and
                                    will soon be archived"</em
                                >
                                (users directed to Flower)
                            </td>
                        </tr>
                        <tr>
                            <td><strong>FeatureCloud</strong></td>
                            <td>
                                <em
                                    >"all applications and the entire
                                    architecture … are open source"</em
                                >
                            </td>
                            <td>
                                no public source located for frontend, backend
                                or controller. Cannot be self-hosted
                            </td>
                        </tr>
                    </tbody>
                </table>
                <p class="small muted" style="margin-top: 0.5em">
                    None of this is a criticism of the work: a property the
                    paper presents as <em>central</em> has been refactored away
                    in two cases, and the citation has not caught up. It is the
                    reason the table is scored from source.
                </p>
            </div>
            <aside class="notes">
                (~1 min 30) This is the payoff of the code-first scoring rule,
                and it's the part of the paper I'd defend hardest, but make sure
                it lands as a statement about method, not about anyone's
                engineering. Substra's published contribution was ledger-backed
                traceability; it was removed in a breaking change and the paper
                is still cited for it. FEDn reset to a client SDK. Both teams
                made deliberate product decisions, and OpenFL's README is
                explicit that it is sunsetting in favour of Flower; they are
                consolidating the ecosystem themselves, which is a good thing
                for users. The lesson is about our instrument: if you score from
                publications you are describing a field that no longer exists.
                Do not editorialise about the platforms; the version and commit
                are in the dossier so anyone can check.
            </aside>
        </section>

        <!-- 22 · Health table — stars measure history, not life.
             Restricted to the peer group from slide 19: the systems a hospital
             weighs against FLIP, not the engines and analytics layers (those
             are consumers or neighbours at another layer, and mixing them in
             was the defect this slide used to have).

             Aggregation rule, stated on the slide because it changes the
             numbers: where a platform spans repositories, STARS are the highest
             single repository and COMMITS are summed across all of them. A
             maximum for stars because a star is a person's judgement about one
             repository and adding two repos' stars would count two different
             audiences as one; a sum for commits because work in any repository
             of the platform is work on the platform.
               · CODA is 32 repositories (its whole GitHub organisation; the
                 study's unit of analysis is the organisation, see slide 4).
               · PHT-meDIC / PADME is the 9 repositories documented on slide 4.
                 The PHT-Medic organisation holds 26 in total, but the other 17
                 are not part of the platform the comparison scored.
               · Every other peer-group row is a single repository.
               · GenoMed4All and INCISIVE publish no repository at all, so their
                 cells are dashes rather than zeros: a zero would read as "we
                 measured it and found nothing", which is false.

             Sources and window: commits and authors from the GitHub API on the
             default branch, 2025-09-06 to 2026-09-06 (the paper's stated
             measurement date), read per repository and then aggregated by the
             rule above. Stars are the GitHub API's current value. Citations are
             OpenAlex for the row's paper (successor paper where a row merges
             two projects, e.g. Kaapana over JIP) — CODA's 21 reproduces the
             figure already in the deck, which is what confirmed the source.

             Two corrections carried here from the local clones: the clones of
             CODA, PHT-MEDIC and FeTS are stale snapshots (HEADs from 2023-2024),
             so a clone-based count reports "0 commits in 12 months" for reasons
             that have nothing to do with the project. Counted from the API, the
             zeros are real but their dates are not "stable since 2025": CODA's
             last commit is 2023-12-22, PHT-meDIC's 2024-01-11, FeTS's
             2023-10-13. The clone-based method is sound where the clone is
             fresh, and reproduces the deck's older figures exactly (Flower
             1,686, FLARE 984, vantage6 1,177). -->
        <section>
            <div class="eyebrow">Project health</div>
            <h2>Health, not popularity</h2>
            <div class="slide-body">
                <p class="small" style="margin: 0 0 0.3em">
                    Stars measure accumulated attention; contributors and commits
                    measure who is still working. Where a platform spans
                    repositories, stars are the
                    <em>highest single repository</em> and commits are
                    <em>summed across all of them</em>. CHAOSS metric families;
                    citations via OpenAlex, Sep 2026.
                </p>
                <table style="font-size: 0.46em">
                    <thead>
                        <tr>
                            <th>Platform</th>
                            <th>Repos</th>
                            <th>Stars</th>
                            <th>Citations</th>
                            <th>Active people</th>
                            <th>Commits 12 mo</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td>FeTS</td>
                            <td>1</td>
                            <td><strong>95</strong></td>
                            <td>384</td>
                            <td>0</td>
                            <td>0</td>
                            <td class="muted">No commits since Oct 2023</td>
                        </tr>
                        <tr>
                            <td>CODA</td>
                            <td>32</td>
                            <td>2</td>
                            <td>21</td>
                            <td>0</td>
                            <td>0</td>
                            <td class="muted">No commits since Dec 2023</td>
                        </tr>
                        <tr>
                            <td>PHT-meDIC / PADME</td>
                            <td>9</td>
                            <td>9</td>
                            <td>11</td>
                            <td>0</td>
                            <td>0</td>
                            <td class="muted">No commits since Jan 2024</td>
                        </tr>
                        <tr>
                            <td>JIP / Kaapana</td>
                            <td>1</td>
                            <td>274</td>
                            <td>7</td>
                            <td>11</td>
                            <td>780</td>
                            <td class="green">Active</td>
                        </tr>
                        <tr>
                            <td>MedPerf</td>
                            <td>1</td>
                            <td>172</td>
                            <td>171</td>
                            <td>4</td>
                            <td>25</td>
                            <td class="green">Active</td>
                        </tr>
                        <tr>
                            <td>NeuroFLAME</td>
                            <td>1</td>
                            <td>3</td>
                            <td>—</td>
                            <td>4</td>
                            <td>238</td>
                            <td class="green">Active</td>
                        </tr>
                        <tr>
                            <td>FLA<sup>3</sup></td>
                            <td>1</td>
                            <td>0</td>
                            <td>—</td>
                            <td>1</td>
                            <td>7</td>
                            <td class="green">New, single author</td>
                        </tr>
                        <tr>
                            <td>GenoMed4All<sup class="accent">⋆</sup></td>
                            <td class="muted">—</td>
                            <td class="muted">—</td>
                            <td>2</td>
                            <td class="muted">—</td>
                            <td class="muted">—</td>
                            <td class="muted">Publication only</td>
                        </tr>
                        <tr>
                            <td>INCISIVE<sup class="accent">⋆</sup></td>
                            <td class="muted">—</td>
                            <td class="muted">—</td>
                            <td>27</td>
                            <td class="muted">—</td>
                            <td class="muted">—</td>
                            <td class="muted">Publication only</td>
                        </tr>
                        <tr class="is-ours">
                            <td><strong>FLIP (ours)</strong></td>
                            <td>1</td>
                            <td>53</td>
                            <td class="muted">—</td>
                            <td>20</td>
                            <td><strong>4,096</strong></td>
                            <td class="green">Active</td>
                        </tr>
                    </tbody>
                </table>
                <p class="small muted" style="margin-top: 0.4em">
                    Commits and authors: default branch, 2025-09-06 →
                    2026-09-06. Stars and citations read Sep 2026.
                </p>
                <p class="small" style="margin-top: 0.45em">
                    <strong>FeTS carries 95 stars and 384 citations and has not
                    been committed to since October 2023</strong>; CODA and
                    PHT-meDIC went quiet the same year while their papers keep
                    being cited. FLIP has 53 stars and more commits than the rest
                    of the table combined. Popularity and momentum are different
                    signals, and a chooser needs both.
                </p>
            </div>
            <aside class="notes">
                (~1 min 30) Same lesson as the popularity trap, arriving from the
                other direction; that is why it lands. Read it as a caution about
                the metric, not about the projects. Three of these are dormant and
                richly cited: FeTS is the sharpest case, 384 citations and 95
                stars with no commit since October 2023, and CODA's paper is
                recent while its last commit is December 2023. That is a
                legitimate state for hospital software — a platform that stops
                changing because it works is a good outcome for the site running
                it — and some of these teams moved to successor projects. What
                the numbers cannot tell you is whether anyone will answer an
                issue next month. Our own row is the mirror image, and don't
                oversell it: 4,096 commits is a team working full time on one
                product, not a virtue in itself, and a young project accumulates
                commits that a mature one has no reason to generate. If someone
                defends a dormant-but-deployed platform, agree: the point is to
                measure both, not to rank them. Two method notes if asked: this
                table is the peer group only, so no engines and no analytics
                layers; and where a platform spans repositories the stars are the
                highest single repository while the commits are summed — CODA is
                32 repositories, PHT-meDIC nine, which is why their zeros carry
                weight rather than being one quiet repo. Citations are OpenAlex
                (September 2026) for the row's paper, cross-checked against
                Crossref on the DOI-bearing ones; treat the column as an
                order-of-magnitude signal, not a ranking. FLIP's citation cell
                stays blank: ours is a preprint under review, and an empty cell is
                more honest than a zero. FLA³'s is blank for the opposite reason,
                a 2026 journal article that is not indexed yet.
            </aside>
        </section>

<!-- 23 · Where FLIP stands. Complementarity first, then our own limit:
         the point is that different platforms serve different needs, and that
         we take specific things from the ones that solved a problem better. -->
        <section>
            <div class="eyebrow">Reading the table honestly</div>
            <h2>Different platforms, different needs</h2>
            <div class="slide-body">
                <div class="cols tight" style="--n: 2">
                    <div class="panel flip">
                        <h3>What FLIP is designed for</h3>
                        <ul class="checklist small">
                            <li>
                                Cohort query in OMOP + MI-CDM, aggregate-only
                                return
                            </li>
                            <li>
                                Resolving to DICOM retrieved on demand from each
                                site's clinical PACS
                            </li>
                            <li>One platform contract over two FL engines</li>
                            <li>
                                Per-project, per-site approval that gates
                                execution
                            </li>
                        </ul>
                        <p class="small muted">
                            A statement about orchestration, not about being
                            better.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Different problems, solved well elsewhere</h3>
                        <ul class="pluslist small">
                            <li>
                                <strong>DataSHIELD, vantage6, MEDIATA</strong>
                                federated analytics without imaging; lighter
                                than FLIP where the study is tabular
                            </li>
                            <li>
                                <strong>FLARE, Flower</strong>: research speed;
                                your own aggregation, no platform to onboard
                            </li>
                            <li>
                                <strong>Swarm Learning, NEBULA</strong>: truly
                                peer-to-peer, where no site will host a hub
                            </li>
                            <li>
                                <strong>Kaapana</strong>: DICOM-native workflow;
                                the best fit when the estate is already imaging
                            </li>
                        </ul>
                    </div>
                </div>
                <div class="flip-note">
                    <span class="flip-tag">Our limit</span> The
                    <em>authority</em>
                    to set approval is one platform-wide permission, so an
                    operator can approve for any site; found by auditing our own
                    code.
                </div>
            </div>
            <aside class="notes">
                (~2 min) The most important slide for this audience, because
                it's the one that decides whether the table reads as evidence or
                as marketing. Lead with the right-hand column and mean it: these
                are not failed versions of FLIP, they are answers to different
                questions, and several are the better answer for questions we
                don't serve. A tabular epidemiology study does not need a PACS
                path; a methods group iterating on aggregation does not want a
                platform to onboard. If a hospital asks us to do something on
                this list, the honest answer is to point at the platform that
                does it well. Then read the bottom box out in full; that is the
                limit we found in our own approval model, and it is in the
                paper.
            </aside>
        </section>

        <!-- 24 · What we take from them. Keeps the comparison from reading as an audit.
        <section>
            <div class="eyebrow">What we take from this</div>
            <h2>Three things worth copying</h2>
            <div class="slide-body">
                <div class="cols" style="--n: 3">
                    <div class="panel">
                        <h3>Enforce policy in the request path</h3>
                        <p class="small">
                            <strong>Kaapana</strong> puts Traefik forwardAuth in
                            front of an Open Policy Agent decision, so a blocked
                            request never reaches the application. FLIP renders
                            the decision in the page; theirs is enforced a layer
                            earlier.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Refuse an unapproved plan by hash</h3>
                        <p class="small">
                            <strong>Fed-BioMed</strong> hashes the training
                            plan, stores it node-side as PENDING/APPROVED, then
                            rejects a mismatch before doing any work, close to
                            what our per-site approval should become.
                        </p>
                    </div>
                    <div class="panel">
                        <h3>Ship a typed audit record</h3>
                        <p class="small">
                            <strong>Flower</strong> emits a
                            <code>LogEntry</code> with actor, event and status;
                            <strong>Substra</strong> and
                            <strong>vantage6</strong> persist lineage per task
                            and run. Worth aligning formats rather than
                            inventing ours alone.
                        </p>
                    </div>
                </div>
                <div class="flip-note">
                    <span class="flip-tag">Ecosystem</span> Three of the four
                    design improvements we would make next came from reading
                    someone else's code.
                </div>
            </div>
            <aside class="notes">
                (~1 min 15) This is the slide that stops the comparison reading
                as "we came first". Each of these is a concrete, adoptable
                improvement with a named source, not flattery. Two of them are
                about enforcement happening earlier than ours; that is the same
                lesson arriving twice, from Kaapana at the network layer and
                from Fed-BioMed at the node, which is a useful thing to say out
                loud. If someone asks why we don't just adopt them wholesale:
                our constraint is the hub-and-spoke model the hospitals signed
                up for, so it is a design trade, but the enforcement point is
                separable and worth taking.
            </aside>
        </section>
 -->
        <!-- 26 · Close. Stays up through Q&A.
             Deliberately just a thank-you. It used to carry three takeaway
             panels and a provenance line; the takeaways are made by the argument
             earlier in the talk, and a closing slide piled with text competes
             with the Q&A it exists to stay up through. The eyebrow repeats the
             title slide's so the talk bookends itself. Nothing was dropped
             silently: the artefact pointers (protocol, queries, inclusion
             criteria, verdict files) and the closing line about the field moved
             into the speaker notes, where they answer a question instead of
             closing a talk. -->
        <section class="title-slide center">
            <div class="eyebrow">Federated Learning Platform Comparison</div>
            <h1>Thank you</h1>
            <p class="subtitle">Questions welcome</p>
            <aside class="notes">
                (~30 s, then Q&amp;A) Leave this up. Give the closing line before
                the questions start: every platform in these tables is the best
                answer to some question, and FLIP's question is cohort query and
                on-demand imaging governed per site — the clearest improvements
                to our own row came from reading the others rather than from
                defending it. If anyone asks for the artefacts, they are released
                with the paper: the search protocol and queries, the inclusion
                criteria, and the per-platform verdict files. The likely
                questions: why not PRISMA (we screened all four pools against
                the five criteria and released the verdict files, but that is
                not a systematic-review-grade literature pass; we do not claim
                it is); why include closed platforms (tier B, labelled,
                excluding them understates what a hospital gets offered);
                whether the engine layer belongs in the same table at all (yes,
                because we consume it, and the governance columns are exactly
                where the engines stop); and which platform they should use
                (depends on the question; that is the point of the slide before
                this one).
            </aside>
        </section>

    </RevealDeck>
</template>

<style>
/* Deck-local: the two capability matrices are 10 columns of mostly single
   glyphs. Auto table layout sizes them from the widest cell and overflows the
   slide, so the widths are declared and the marker columns centred by position
   rather than by a class on ~180 cells. Selectors are keyed to .cap-table so
   nothing here reaches the other tables in the deck. */
.reveal.deck-theme .cap-table {
    table-layout: fixed;
}
.reveal.deck-theme .cap-table th,
.reveal.deck-theme .cap-table td {
    padding: 0.26em 0.35em;
}
.reveal.deck-theme .cap-table th:nth-child(n + 3),
.reveal.deck-theme .cap-table td:nth-child(n + 3) {
    text-align: center;
}
/* .cov carries a right margin for inline use; centred in a cell that reads as
   off-centre. */
.reveal.deck-theme .cap-table .cov {
    margin-right: 0;
}
/* Deck-local: on the complementarity slide the two lists ARE the content, so
   the per-item gap is tightened to keep the closing note on the canvas. */
.reveal.deck-theme .cols.tight li {
    margin-bottom: 0.15em;
}
/* Our own row, in a table of other people's work. */
.reveal.deck-theme table tr.is-ours td {
    background: var(--flip-col-bg);
}
/* Platform-name links in the capability tables. Reveal's default link colour is
   too loud for a 15-row table read at a glance: the eye should still land on the
   coverage glyphs, not on a column of blue. These inherit the cell's colour and
   carry a faint underline instead, so the affordance is visible up close (and on
   hover) without competing with the data. */
.reveal.deck-theme .cap-table td a {
    color: inherit;
    text-decoration: none;
    border-bottom: 1px dotted currentColor;
    opacity: 0.92;
}
.reveal.deck-theme .cap-table td a:hover,
.reveal.deck-theme .cap-table td a:focus-visible {
    color: var(--accent-cyan, var(--accent));
    border-bottom-style: solid;
    opacity: 1;
}
/* Deck-local: a neutral list marker. .crosslist puts a pink ✗ on every line,
   which is the wrong signal when the line is describing another platform's
   strength. The comparison names trade-offs, it does not strike them out. */
.reveal.deck-theme .pluslist {
    list-style: none;
    margin-left: 0;
}
.reveal.deck-theme .pluslist li::before {
    content: "\2192";
    color: var(--accent-cyan, var(--accent));
    font-weight: 700;
    margin-right: 0.5em;
}
</style>
