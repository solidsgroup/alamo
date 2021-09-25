/* IMPORTANT
 * This snapshot file is auto-generated, but designed for humans.
 * It should be checked into source control and tracked carefully.
 * Re-generate by setting TAP_SNAPSHOT=1 and running tests.
 * Make sure to inspect the output below.  Do not ignore changes!
 */
'use strict'
exports[`test/cmd.js TAP -h --help prints usage > --help output 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "\\nusage: mkdirp [DIR1,DIR2..] {OPTIONS}\\n\\n  Create each supplied directory including any necessary parent directories\\n  that don't yet exist.\\n\\n  If the directory already exists, do nothing.\\n\\nOPTIONS are:\\n\\n  -m<mode>       If a directory needs to be created, set the mode as an octal\\n  --mode=<mode>  permission string.\\n\\n  -v --version   Print the mkdirp version number\\n\\n  -h --help      Print this helpful banner\\n\\n  -p --print     Print the first directories created for each path provided\\n\\n  --manual       Use manual implementation, even if native is available\\n\\n",
}
`

exports[`test/cmd.js TAP -v --version prints version > --version output 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "4.2.0-69.lol\\n",
}
`

exports[`test/cmd.js TAP failures > expect resolving Promise 1`] = `
Array [
  Object {
    "code": 1,
    "signal": null,
    "stderr": "nope\\n",
    "stdout": "",
  },
  Object {
    "code": 1,
    "signal": null,
    "stderr": "fail\\n  code: EFAIL\\n",
    "stdout": "",
  },
]
`

exports[`test/cmd.js TAP invalid mode > expect resolving Promise 1`] = `
Object {
  "code": 1,
  "signal": null,
  "stderr": "invalid mode argument: --mode=XYZ\\nMust be an octal number.\\n",
  "stdout": "",
}
`

exports[`test/cmd.js TAP make dir named --help > expect resolving Promise 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "--help 0\\n",
}
`

exports[`test/cmd.js TAP making dirs > expect resolving Promise 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "",
}
`

exports[`test/cmd.js TAP manual > expect resolving Promise 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "MANUAL a 0\\nMANUAL b/c/d 0\\n",
}
`

exports[`test/cmd.js TAP no dirs -> stderr usage > expect resolving Promise 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "\\nusage: mkdirp [DIR1,DIR2..] {OPTIONS}\\n\\n  Create each supplied directory including any necessary parent directories\\n  that don't yet exist.\\n\\n  If the directory already exists, do nothing.\\n\\nOPTIONS are:\\n\\n  -m<mode>       If a directory needs to be created, set the mode as an octal\\n  --mode=<mode>  permission string.\\n\\n  -v --version   Print the mkdirp version number\\n\\n  -h --help      Print this helpful banner\\n\\n  -p --print     Print the first directories created for each path provided\\n\\n  --manual       Use manual implementation, even if native is available\\n\\n",
  "stdout": "",
}
`

exports[`test/cmd.js TAP noisily > expect resolving Promise 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "a 0\\nb/c/d 0\\n",
}
`

exports[`test/cmd.js TAP print modes > expect resolving Promise 1`] = `
Object {
  "code": 0,
  "signal": null,
  "stderr": "",
  "stdout": "a 509\\n",
}
`
