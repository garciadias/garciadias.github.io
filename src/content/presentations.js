// Registry of slide decks. Each entry feeds the /presentations gallery and the
// fullscreen /presentations/:id deck route. `deck` is a lazy import of the Vue
// component holding the reveal.js <section> slides.
//
// To add a new presentation: create src/decks/<Name>Deck.vue, then add an entry
// here. Nothing else needs editing.
//
// Set `unlisted: true` to publish a deck without showing it on the public
// /presentations gallery. Unlisted decks still work at their direct
// /presentations/:id URL and appear on the secret /presentations/list index.
export const presentations = [
  {
    id: 'iaa-so-chemical-tagging-2026',
    title: 'Chemical tagging: finding lost star clusters',
    subtitle:
      'Unsupervised learning on stellar spectra, from K-means to EVoC · IAA-SO School on AI/ML in Astronomy 2026',
    date: '2026',
    venue: 'IAA-SO School on AI/ML in Astronomy 2026 · Unsupervised Learning pillar',
    cover: `${import.meta.env.BASE_URL}presentations/iaa-so-chemical-tagging-2026/benchmark_grid.png`,
    description:
      'A ~90-minute lecture plus hands-on for the Unsupervised Learning pillar of the IAA-SO school. ' +
      'Stars born together share a chemical fingerprint; clusters dissolve but chemistry does not, so ' +
      'can we reconstruct them from the spectra alone? Builds the clustering toolkit in order — ' +
      'K-means, KNN, DBSCAN, HDBSCAN*, PLSCAN, t-SNE, UMAP, EVoC — and applies it to 25 open and ' +
      'globular clusters in APOGEE DR19 + Gaia DR3, reproducing the opposite verdicts of Kos et al. ' +
      '2017 and Garcia-Dias et al. 2019 before extending them. The result: a masked spectral ' +
      'autoencoder trained with no labels at all separates clusters better than the ASPCAP ' +
      'abundances (homogeneity 0.79/0.87 vs 0.48 on a matched sample), and the win survives ' +
      'removing the globular that dominates the sample. Includes the batch-effect and ' +
      'seed-stability controls behind those numbers, and a student assignment: take a cluster, ' +
      'beat the baseline.',
    tags: [
      'Unsupervised Learning',
      'Astronomy',
      'Chemical Tagging',
      'Clustering',
      'Self-Supervised Learning',
      'APOGEE',
      'IAA-SO 2026'
    ],
    deck: () => import('@/decks/IaaSoChemicalTaggingDeck.vue')
  },
  {
    id: 'flare-day-2026',
    title: 'FLIP: an open-source federated learning platform for healthcare',
    subtitle: 'From multi-institutional research to real NHS deployment · NVIDIA FLARE Day 2026',
    date: '16 September 2026',
    venue: 'NVIDIA FLARE Day 2026 · US + EMEA main event',
    cover: `${import.meta.env.BASE_URL}presentations/flare-day-2026/intro-globe-poster.jpg`,
    description:
      'A 20-minute talk for NVIDIA FLARE Day on what sits between a federated learning framework ' +
      'and a hospital. FLARE solved the federated core; every collaboration still rebuilds cohort ' +
      'definition, imaging retrieval, per-site approval, scheduling and audit by hand. FLIP — an ' +
      'open-source project by the London AI Centre, King\'s College London and Guy\'s and St Thomas\' ' +
      'NHS Foundation Trust — is that layer, solved once. Covers where FLIP sits on the NVIDIA stack ' +
      '(what it delegates to FLARE, and MONAI in both directions: unmodified bundles training as a ' +
      '3D segmentation job type, trained models leaving as MONAI Application Packages), the single ' +
      'outbound HTTPS connection it asks an NHS network for, cohort query and per-site approval ' +
      'shown in the product, the six sovereignty guarantees and compliance frameworks behind it, ' +
      'a live UK ⇄ Thailand run and a standing clinical study across three trusts, the AWS Landing ' +
      'Zone hub, and three lessons that cost us months: certificate rotation, coding heterogeneity ' +
      'ahead of statistical heterogeneity, and governance as the long pole. Closes with a short cut ' +
      'of our 30-platform capability audit, dotted against the FLARE Day programme.',
    tags: ['Federated Learning', 'FLIP', 'NVIDIA FLARE', 'MONAI', 'NHS', 'Open Source', 'FLARE Day'],
    deck: () => import('@/decks/FlareDay2026Deck.vue')
  },
  {
    id: 'flip-platform-comparison-amigo',
    title: 'Comparing FL platforms — fairly',
    subtitle: 'Inclusion criteria, a reproducible search, and a 30-platform capability audit scored on code, not claims',
    date: 'September 2026',
    venue: 'AMIGO team meeting',
    cover: `${import.meta.env.BASE_URL}presentations/flip-maturity-pitch-2026/flip-architecture-symmetric.png`,
    description:
      'The methods talk behind the comparison tables in our MICCAI 2026 / DeCaF paper. How do you ' +
      'decide which federated learning platforms belong in a comparison, and how do you score them ' +
      'so a reviewer with repository access cannot take the table apart? Covers the eligibility ' +
      'criteria and evidence tiers, a reproducible four-source search over PubMed, GitHub, arXiv ' +
      'and medRxiv, and the measured recall of each — PubMed recovers 13 of 30 and misses six of ' +
      'the eight FL engines outright; GitHub reaches 20 of the 22 with a public repository, but ' +
      'only once the metadata-poor platforms are queried by name; together the four recover all ' +
      '30. Then the capability audit itself across 30 platforms and nine columns, every cell ' +
      're-read from a local clone at a recorded commit rather than from publications — which ' +
      'moved 48 cells and is how we found that two platforms have lost the headline property ' +
      'their papers are still cited for. Includes where FLIP is beaten, and the limit we found ' +
      'in our own approval model.',
    tags: ['Federated Learning', 'Platform Comparison', 'Systematic Search', 'Research Methods', 'FLIP', 'MICCAI 2026'],
    deck: () => import('@/decks/PlatformComparisonDeck.vue')
  },
  {
    id: 'fla3-governance-federated-learning',
    title: 'What can FLIP learn from FLA³?',
    subtitle: 'A FLIP-eyed read of FLA³ — Federated Learning with Authentication · Authorisation · Accounting',
    date: '2 Jun 2026',
    venue: 'AMIGO team meeting',
    cover: `${import.meta.env.BASE_URL}presentations/fla3/x4.png`,
    description:
      'A walkthrough of the FLA³ paper (arXiv 2603.10063) through one question: what can our FLIP platform take from it? FLA³ enforces Authentication, Authorisation & Accounting as a runtime control plane — per-round XACML policy decisions, fail-closed, with cryptographically signed audit. Its thesis: in regulated healthcare the breach is unauthorised computation, not data movement — and governance costs nothing in accuracy (federation matches centralised on the INTERVAL iron-deficiency task, and lifts the weakest sites most). Every slide carries an honest FLIP verdict: where we already match, and the upgrades worth weighing.',
    tags: ['Federated Learning', 'Governance', 'Healthcare AI', 'FLIP', 'FLA³'],
    deck: () => import('@/decks/Fla3Deck.vue')
  },
  {
    id: 'miccai_decaf_2026_draft',
    title: 'Federated chest X-ray learning, UK ⇄ Thailand',
    subtitle: 'Cross-continental federated fine-tuning of a CXR foundation model with FLIP — MICCAI 2026 / DECAF working draft',
    date: '2026',
    venue: 'MICCAI 2026 · DECAF (draft)',
    description:
      'A working-draft walkthrough of a UK–Thailand proof of concept: deploying FLIP — our Federated Learning & Interoperability Platform, validated at NHS scale — in a setting it was not designed for, a live cross-continental, cross-sector federation linking UK academia with a private Thai hospital group. Not a new algorithm: a reproducible federated fine-tuning workflow for a pretrained chest X-ray foundation model, a four-arm comparison (zero-shot · local-only · federated · centralised) over locally generated synthetic data, and an evaluation that pairs predictive metrics with system- and deployment-level measurements. Unlisted draft for internal review.',
    tags: ['Federated Learning', 'Chest X-ray', 'Foundation Models', 'FLIP', 'MICCAI 2026', 'Draft'],
    unlisted: true,
    deck: () => import('@/decks/MiccaiDecaf2026Deck.vue')
  },
  {
    id: 'pydata-london-flip-lightning-2026',
    title: 'AI on the NHS: can federated learning preserve patient privacy?',
    subtitle: 'A 5-minute case for sending the model to the data, not the other way round',
    date: '2026',
    venue: 'PyData London · Lightning talk',
    cover: `${import.meta.env.BASE_URL}presentations/flip-maturity-pitch-2026/flip-architecture-symmetric.png`,
    description:
      'Federated learning inverts the logic of big-data lakes for training AI models: instead of bringing ' +
      'data to models, we send the models to be trained where the data is. A fast, informal cut of FLIP — ' +
      'the platform the AI Centre for Value-Based Healthcare runs in production across the NHS and a ' +
      'Thai hospital group — built for a PyData lightning-talk slot: the problem, the one idea, the ' +
      'architecture in one diagram, real production numbers, and the open-source Python stack underneath it.',
    tags: ['Federated Learning', 'Privacy', 'NHS', 'PyData', 'Lightning Talk', 'FLIP'],
    deck: () => import('@/decks/PydataLondonFlipLightningDeck.vue')
  },
  {
    id: 'flip-reusable-fl-platform-2026',
    title: 'FLIP: a reusable, general-purpose open-source FL platform',
    subtitle: 'From Zero to Hero: Federated AI in Healthcare Systems · St John\'s College, Cambridge',
    date: '2026',
    venue: 'Federated AI in Healthcare Systems · Cambridge',
    cover: `${import.meta.env.BASE_URL}presentations/flip-reusable-fl-platform/intro-globe-poster.jpg`,
    description:
      'A 25-minute workshop talk on what it takes to turn federated learning from a research project ' +
      'into standing infrastructure. The argument: the FL frameworks are ready — Flower and NVIDIA FLARE ' +
      'made the federated core something you can build on — but every collaboration still rebuilds ' +
      'cohort definition, harmonisation, per-site approval, scheduling and audit by hand. FLIP is that ' +
      'layer between the hospital and the framework, solved once, open source and Apache 2.0. Covers the ' +
      'services it decomposes into, the single outbound connection it asks a hospital network for, and ' +
      'where the guardrails execute — then a demo of a real UK ⇄ Thailand deployment: cohort query, ' +
      'per-site approval, and a live federated training run.',
    tags: ['Federated Learning', 'FLIP', 'NHS', 'Open Source', 'Platform', 'Cambridge'],
    deck: () => import('@/decks/FlipReusableFlPlatformDeck.vue')
  }
]

// Decks shown on the public gallery (everything not flagged `unlisted`).
export const listedPresentations = presentations.filter((p) => !p.unlisted)

export function getPresentation(id) {
  return presentations.find((p) => p.id === id)
}
