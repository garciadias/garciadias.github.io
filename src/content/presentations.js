// Registry of slide decks. Each entry feeds the /presentations gallery and the
// fullscreen /presentations/:id deck route. `deck` is a lazy import of the Vue
// component holding the reveal.js <section> slides.
//
// To add a new presentation: create src/decks/<Name>Deck.vue, then add an entry
// here. Nothing else needs editing.
export const presentations = [
  {
    id: 'fla3-governance-federated-learning',
    title: 'What can FLIP learn from FLA³?',
    subtitle: 'A FLIP-eyed read of FLA³ — Federated Learning with Authentication · Authorisation · Accounting',
    date: '2 Jun 2026',
    venue: 'AMIGO team meeting',
    cover: `${import.meta.env.BASE_URL}presentations/fla3/cover.png`,
    description:
      'A walkthrough of the FLA³ paper (arXiv 2603.10063) through one question: what can our FLIP platform take from it? FLA³ enforces Authentication, Authorisation & Accounting as a runtime control plane — per-round XACML policy decisions, fail-closed, with cryptographically signed audit. Its thesis: in regulated healthcare the breach is unauthorised computation, not data movement — and governance costs nothing in accuracy (federation matches centralised on the INTERVAL iron-deficiency task, and lifts the weakest sites most). Every slide carries an honest FLIP verdict: where we already match, and the upgrades worth weighing.',
    tags: ['Federated Learning', 'Governance / AAA', 'XACML PDP', 'Healthcare AI', 'FLIP'],
    deck: () => import('@/decks/Fla3Deck.vue')
  }
]

export function getPresentation(id) {
  return presentations.find((p) => p.id === id)
}
