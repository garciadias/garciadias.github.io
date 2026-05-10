import { createRouter, createWebHashHistory } from 'vue-router'

const routes = [
  {
    path: '/',
    name: 'home',
    component: () => import('@/views/HomeView.vue'),
    meta: { title: 'Home' }
  },
  {
    path: '/cv',
    name: 'cv',
    component: () => import('@/views/CVView.vue'),
    meta: { title: 'Career' }
  },
  {
    path: '/projects',
    name: 'projects',
    component: () => import('@/views/ProjectsView.vue'),
    meta: { title: 'Projects' }
  },
  {
    path: '/publications',
    name: 'publications',
    component: () => import('@/views/PublicationsView.vue'),
    meta: { title: 'Publications' }
  },
  {
    path: '/talks',
    name: 'talks',
    component: () => import('@/views/TalksView.vue'),
    meta: { title: 'Talks & Videos' }
  },
  {
    path: '/:pathMatch(.*)*',
    name: 'not-found',
    component: () => import('@/views/NotFoundView.vue'),
    meta: { title: 'Not found' }
  }
]

const router = createRouter({
  history: createWebHashHistory(),
  routes,
  scrollBehavior(to, from, savedPosition) {
    if (savedPosition) return savedPosition
    return { top: 0 }
  }
})

router.afterEach((to) => {
  const base = 'Rafael Garcia-Dias'
  document.title = to.meta?.title ? `${to.meta.title} — ${base}` : base
})

export default router
