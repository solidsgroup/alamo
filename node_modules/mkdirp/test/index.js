const t = require('tap')
const mkdirp = require('../')

t.test('module shape', t => {
  t.isa(mkdirp, Function)
  t.isa(mkdirp.sync, Function)
  t.isa(mkdirp.manual, Function)
  t.isa(mkdirp.manualSync, Function)
  t.isa(mkdirp.native, Function)
  t.isa(mkdirp.nativeSync, Function)
  t.end()
})

t.test('basic making of dirs should work', t => {
  const dir = t.testdir({ a: {} })
  const {statSync, mkdir, mkdirSync} = require('fs')
  const check = d => t.ok(statSync(d).isDirectory())
  t.equal(mkdirp.sync(`${dir}/a/sync`), `${dir}/a/sync`)
  check(`${dir}/a/sync`)
  t.equal(mkdirp.sync(`${dir}/a/sync`), undefined)

  t.equal(mkdirp.manualSync(`${dir}/a/manual-sync`), `${dir}/a/manual-sync`)
  check(`${dir}/a/manual-sync`)
  t.equal(mkdirp.manualSync(`${dir}/a/manual-sync`), undefined)

  t.equal(mkdirp.nativeSync(`${dir}/a/native-sync`), `${dir}/a/native-sync`)
  check(`${dir}/a/native-sync`)
  t.equal(mkdirp.nativeSync(`${dir}/a/native-sync`), undefined)

  // override to force the manual option
  const myMkdir = (path, opts, cb) => mkdir(path, opts, cb)
  const myMkdirSync = (path, opts) => mkdirSync(path, opts)
  const opts = { mkdir: myMkdir, mkdirSync: myMkdirSync }
  t.equal(mkdirp.sync(`${dir}/a/custom-sync`, opts), `${dir}/a/custom-sync`)
  check(`${dir}/a/custom-sync`)
  t.equal(mkdirp.sync(`${dir}/a/custom-sync`, opts), undefined)

  return Promise.all([
    mkdirp(`${dir}/a/async`),
    mkdirp.manual(`${dir}/a/manual-async`),
    mkdirp.native(`${dir}/a/native-async`),
    mkdirp(`${dir}/a/custom-async`, opts),
  ]).then(made => {
    t.strictSame(made, [
      `${dir}/a/async`,
      `${dir}/a/manual-async`,
      `${dir}/a/native-async`,
      `${dir}/a/custom-async`,
    ])
    check(`${dir}/a/async`)
    check(`${dir}/a/manual-async`)
    check(`${dir}/a/native-async`)
    check(`${dir}/a/custom-async`)
    return Promise.all([
      mkdirp(`${dir}/a/async`),
      mkdirp.manual(`${dir}/a/manual-async`),
      mkdirp.native(`${dir}/a/native-async`),
      mkdirp(`${dir}/a/custom-async`, opts),
    ])
  }).then(made => t.strictSame(made, [undefined, undefined, undefined, undefined]))
})
