const t = require('tap')
const requireInject = require('require-inject')

const {basename, posix} = require('path')
const {promisify} = require('util')
const fs = require('fs')

const statAsync = (path) =>
  basename(path) === 'error'
    ? Promise.reject(new Error('not a real error'))
    : promisify(fs.stat)(path)

const statSync = path => {
  if (basename(path) === 'error')
    throw new Error('not a real error')
  else
    return fs.statSync(path)
}

const {findMade, findMadeSync} = requireInject('../lib/find-made.js', {
  path: posix,
})

t.test('find what dir will be made', t => {
  const dir = t.testdir({
    file: 'txt',
    subdir: {},
  })

  const o = {statAsync, statSync}

  t.equal(findMadeSync(o, `${dir}/subdir/x/y/z`), `${dir}/subdir/x`)
  t.equal(findMadeSync(o, `${dir}/subdir`), undefined)
  t.equal(findMadeSync(o, `${dir}/file/x/y/z`), undefined)
  t.equal(findMadeSync(o, `${dir}/file`, `${dir}/file/x`), undefined)
  t.equal(findMadeSync(o, `${dir}/subdir/error`), undefined)
  t.equal(findMadeSync(o, '/', '/'), undefined)
  return Promise.all([
    findMade(o, `${dir}/subdir/x/y/z`),
    findMade(o, `${dir}/subdir`),
    findMade(o, `${dir}/file/x/y/z`),
    findMade(o, `${dir}/file`, `${dir}/file/x`),
    findMade(o, `${dir}/subdir/error`),
    findMade(o, '/', '/'),
  ]).then(made => t.strictSame(made, [
    `${dir}/subdir/x`,
    undefined,
    undefined,
    undefined,
    undefined,
    undefined,
  ]))
})
