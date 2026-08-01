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
      'architecture in one diagram, real production numbers, and the open-source Python stack underneath ' +
      'it. Draft — pending final gifs and confirmed talk date.',
    tags: ['Federated Learning', 'Privacy', 'NHS', 'PyData', 'Lightning Talk', 'FLIP'],
    unlisted: true,
    deck: () => import('@/decks/PydataLondonFlipLightningDeck.vue')
  }
]

// Decks shown on the public gallery (everything not flagged `unlisted`).
export const listedPresentations = presentations.filter((p) => !p.unlisted)

export function getPresentation(id) {
  return presentations.find((p) => p.id === id)
}
