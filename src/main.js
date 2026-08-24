import { createApp } from 'vue'
import './assets/main.css'

// Two ways in. Normally this is the router-driven site. With ?deck=<id> it is
// present mode: a single deck, no router in the page, so reveal.js owns
// location.hash — which is what makes its speaker view (S) work. See
// views/PresentView.vue for why that matters. The router module is not even
// loaded in present mode, so nothing else can touch the hash.
const presentId = new URLSearchParams(window.location.search).get('deck')

if (presentId) {
  import('./views/PresentView.vue').then(({ default: PresentView }) => {
    createApp(PresentView, { id: presentId }).mount('#app')
  })
} else {
  Promise.all([import('./App.vue'), import('./router')]).then(([{ default: App }, { default: router }]) => {
    const app = createApp(App)
    app.use(router)
    router.isReady().then(() => app.mount('#app'))
  })
}
