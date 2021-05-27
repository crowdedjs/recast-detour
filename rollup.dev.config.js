import babel from '@rollup/plugin-babel';

export default {
  input: './src/Main.js',
  output: {
    file: 'crowded.dev.js',
    format: 'umd',
    name:'crowded'
  },
  plugins: [
    
  ]
};