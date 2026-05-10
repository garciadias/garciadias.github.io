# garciadias.github.io

Personal website of **Dr Rafael Garcia-Dias** — Senior AI Engineer working on federated learning, foundation models and healthcare AI.

Built with **Vue 3 + Vite + Tailwind CSS**, deployed as a static site to **GitHub Pages** via GitHub Actions.

## Local development

```bash
npm install
npm run dev      # local dev server with HMR
npm run build    # static build into ./dist
npm run preview  # preview the production build locally
```

## Project structure

```
public/                 # served as-is at the site root
  static/img/...        # photos, PDFs, publication assets (legacy paths preserved)
src/
  assets/main.css       # global styles + Tailwind layers
  components/           # SiteHeader, SiteFooter, SocialIcon
  composables/          # useDarkMode (light/dark toggle, persisted)
  content/              # all editable content
    profile.js          # bio, headline, metrics, stack, socials
    experience.js       # CV roles + education
    projects.json       # portfolio cards
    publications.json   # selected publications
    talks.json          # talks & embedded YouTube videos
  router/               # Vue Router (hash mode for GitHub Pages)
  views/                # one .vue per page
  App.vue
  main.js
index.html
vite.config.js
tailwind.config.js
.github/workflows/deploy.yml
```

## Adding content

Everything is data-driven — no rebuilding components needed.

- **A new project** → append an entry to `src/content/projects.json`.
- **A new publication** → append to `src/content/publications.json`.
- **A new talk / video** → append to `src/content/talks.json` with a `youtubeId` (the part after `?v=`).
- **CV / experience** → edit `src/content/experience.js`.
- **Profile, headline, stack, socials** → edit `src/content/profile.js`.

Commit and push to `master`/`main` — the GitHub Actions workflow builds and deploys automatically.

## Adding a new section / tab

1. Create `src/views/MySection.vue`.
2. Register it in `src/router/index.js`.
3. Add the route to the `links` array in `src/components/SiteHeader.vue`.

## Deployment

`.github/workflows/deploy.yml` runs on every push to `master` (or `main`):

1. Installs deps with `npm ci`
2. Builds with `npm run build`
3. Publishes `./dist` to GitHub Pages via `actions/deploy-pages@v4`

To enable: in the repo on GitHub go to **Settings → Pages → Build and deployment → Source = GitHub Actions**.

The site uses **hash-mode routing** (`/#/cv`, `/#/projects` …) so it works on GitHub Pages with no 404 fallback hacks.
