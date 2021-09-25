const t = require('tap')
const optsArg = require('../lib/opts-arg.js')
const umask = process.umask()
const mode = 0o777 & (~umask)
const fs = require('fs')

const defFs = {
  fs,
  mkdir: fs.mkdir,
  mkdirSync: fs.mkdirSync,
  stat: fs.stat,
  statSync: fs.statSync,
}

const stat = () => {}
stat.fake = true
const statSync = () => {}
statSync.fake = true

// arg, expect
const cases = {
  null: [null, { mode, ...defFs }],
  false: [false, { mode, ...defFs }],
  undefined: [undefined, { mode, ...defFs }],
  'empty object': [{}, { mode, ...defFs }],
  'numeric mode': [0o775, {mode: 0o775, ...defFs}],
  'string mode': ['775', {mode: 0o775, ...defFs}],
  'empty custom fs': [{fs: {}}, {mode, ...defFs, fs: {}}],
  'custom stat/statSync': [{stat, statSync}, {mode, ...defFs, stat, statSync}],
  'custom fs with stat/statSync': [{fs: {stat, statSync}}, {mode, ...defFs, fs: {stat, statSync}, stat, statSync}],
}

for (const [name, c] of Object.entries(cases)) {
  const [arg, expect] = c
  t.match(optsArg(arg), expect, name)
}

t.throws(() => optsArg(() => {}), TypeError('invalid options argument'))
