const t = require('tap')
const requireInject = require('require-inject')
const {promisify} = require('util')
const {stat, statSync} = require('fs')
const statAsync = promisify(stat)

const {mkdirpNative, mkdirpNativeSync} = requireInject('../lib/mkdirp-native', {
  // just return an indicator that it was called
  '../lib/mkdirp-manual.js': {
    mkdirpManual: () => 'mkdirpManual',
    mkdirpManualSync: () => 'mkdirpManualSync',
  },
  path: require('path').posix
})

const {resolve} = require('path').posix

t.test('mkdirpNative / just calls implementation', t => {
  const opt = {
    mkdirAsync: () => 'mkdirAsync impl',
    mkdirSync: () => 'mkdirSync impl',
  }
  t.equal(mkdirpNative('/', opt), 'mkdirAsync impl')
  t.equal(opt.recursive, true)
  delete opt.recursive
  t.equal(mkdirpNativeSync('/', opt), 'mkdirSync impl')
  t.equal(opt.recursive, true)
  t.end()
})

t.test('mkdirpNative calls impl and returns findMade', t => {
  const opt = {
    mkdirAsync: () => Promise.resolve(),
    mkdirSync: () => undefined,
    statAsync,
    statSync,
  }

  const dir = t.testdir()
  t.equal(mkdirpNativeSync(`${dir}/sync/a/b/c`, opt), `${dir}/sync`)
  return mkdirpNative(`${dir}/async/a/b/c`, opt).then(made =>
    t.equal(made, `${dir}/async`))
})

t.test('ENOENT error falls back to manual', t => {
  const opt = {
    mkdirAsync: () => Promise.reject(Object.assign(new Error('poo'), { code: 'ENOENT' })),
    mkdirSync: () => {
      throw Object.assign(new Error('poo'), { code: 'ENOENT' })
    },
    statAsync,
    statSync,
  }

  const dir = t.testdir()
  t.equal(mkdirpNativeSync(`${dir}/sync/a/b/c`, opt), 'mkdirpManualSync')
  return mkdirpNative(`${dir}/async/a/b/c`, opt).then(made =>
    t.equal(made, 'mkdirpManual'))
})

t.test('other errors are raised to caller', t => {
  const opt = {
    mkdirAsync: () => Promise.reject(Object.assign(new Error('poo'), { code: 'blorg' })),
    mkdirSync: () => {
      throw Object.assign(new Error('poo'), { code: 'blorg' })
    },
    statAsync,
    statSync,
  }

  t.throws(() => mkdirpNativeSync('anything', opt), {code: 'blorg'})
  return t.rejects(mkdirpNative('at/all', opt), {code: 'blorg'})
})
