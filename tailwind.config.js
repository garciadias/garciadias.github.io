/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './src/**/*.{vue,js,ts,jsx,tsx}'],
  darkMode: 'class',
  theme: {
    extend: {
      fontFamily: {
        sans: ['Inter', 'ui-sans-serif', 'system-ui', 'sans-serif'],
        mono: ['JetBrainsMono', 'ui-monospace', 'SFMono-Regular', 'monospace'],
        display: ['"Bai Jamjuree"', 'Inter', 'sans-serif']
      },
      colors: {
        primary: {
          50: '#f0f7ff',
          100: '#e0eefe',
          200: '#bbddfd',
          300: '#7cc1fb',
          400: '#36a3f6',
          500: '#0c87e7',
          600: '#006bc6',
          700: '#0156a0',
          800: '#064a85',
          900: '#0b3e6f',
          950: '#072749'
        }
      },
      typography: {
        DEFAULT: {
          css: {
            maxWidth: 'none'
          }
        }
      }
    }
  },
  plugins: []
}
