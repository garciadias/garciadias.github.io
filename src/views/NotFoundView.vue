<template>
  <div class="container-prose py-24 text-center">
    <p class="font-mono text-sm text-primary-600 dark:text-primary-400">404</p>
    <h1 class="mt-3 text-3xl font-bold text-gray-900 dark:text-white">Page not found</h1>
    <p class="mt-3 text-gray-600 dark:text-gray-400">The page you’re looking for doesn’t exist.</p>
    <router-link
      :to="{ name: 'home' }"
      class="mt-6 inline-flex items-center gap-2 rounded-lg bg-primary-600 hover:bg-primary-700 text-white px-4 py-2 text-sm font-semibold"
      >Back home</router-link
    >
  </div>
</template>
