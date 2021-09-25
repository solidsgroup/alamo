const {basename} = require('path')
const map = base =>
  base === 'index.js' ? 'index.js'
  : base === 'cmd.js' ? 'bin/cmd.js'
  : `lib/${base}`
module.exports = test => map(basename(test))
