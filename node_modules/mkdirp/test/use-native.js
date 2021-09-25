const t = require('tap')
const {mkdir, mkdirSync} = require('fs')
const {useNative, useNativeSync} = require('../lib/use-native.js')

if (!process.env.__TESTING_MKDIRP_NODE_VERSION__) {
  t.spawn(process.execPath, [__filename], {
    env: {
      ...process.env,
      __TESTING_MKDIRP_NODE_VERSION__: 'v10.11.12',
    },
  })

  t.spawn(process.execPath, [__filename], {
    env: {
      ...process.env,
      __TESTING_MKDIRP_NODE_VERSION__: 'v8.9.10',
    },
  })

  // this one has the native impl
  t.equal(useNative({mkdir}), true)
  t.equal(useNative({mkdir: 1243}), false)
  t.equal(useNativeSync({mkdirSync}), true)
  t.equal(useNativeSync({mkdirSync: 1243}), false)
} else {
  t.equal(useNative({mkdir}), false)
  t.equal(useNative({mkdir: 1243}), false)
  t.equal(useNativeSync({mkdirSync}), false)
  t.equal(useNativeSync({mkdirSync: 1243}), false)
}
