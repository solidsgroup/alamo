# One-off changes

Classification: chamber-gpu corpus ledger. Entries are retained as evidence and
must not be treated as a new port's closure or worklist.

Policy: A port worker touching a file with a [NUM] entry must surface it to the user and never apply it silently. Every entry is file-specific GPU-port residue, not a reusable transform.

- src/BC/Operator/Elastic/Expression.H:6-11 | Required GPU-port change to compile or execute BASE functionality on device. | commit fe3844f31579972b174664616e34a0d619c896e6 | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/Constant.H:7-12 | Required GPU-port change to compile or execute BASE functionality on device. | commit f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/Expression.H:34-39 | Required GPU-port change to compile or execute BASE functionality on device. | commit fe3844f31579972b174664616e34a0d619c896e6 | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/Expression.H:44-55 | Required GPU-port change to compile or execute BASE functionality on device. | commit fe3844f31579972b174664616e34a0d619c896e6 | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/Laminate.H:5-10 | Required GPU-port change to compile or execute BASE functionality on device. | commit 98f8ceef437b01144385382ece557d5809e0507f | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/PNG.H:4-9 | Required GPU-port change to compile or execute BASE functionality on device. | commit 1b1b0c1efe224c2a381992d20a5a3b2f18afada2 | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/PSRead.H:6-11 | Required GPU-port change to compile or execute BASE functionality on device. | commit bb4afc9efeb774df2af978989768235511ac6ae2 | file-specific BASE GPU-port residue with no repeatable transform
- src/IC/Trig.H:9-14 | Required GPU-port change to compile or execute BASE functionality on device. | commit fe3844f31579972b174664616e34a0d619c896e6 | file-specific BASE GPU-port residue with no repeatable transform
- src/Integrator/Flame.cpp:595-693#3 | Required GPU-port change to compile or execute BASE functionality on device. | commit adce42a4e9e105f450bbddd4298805d1febdd346 | file-specific BASE GPU-port residue with no repeatable transform
- src/Integrator/PFC.cpp:10-15 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit 9416d16ca695f632d243dc60f153dea57b36752f | file-specific BASE GPU-port residue with no repeatable transform
- src/Model/Propellant/Homogenize.H:8-14 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit c1ec304e2cb201466f7cab623206e7c2e7e9458b | file-specific BASE GPU-port residue with no repeatable transform
- src/Model/Solid/Finite/NeoHookean.H:14-34#2 [NUM] | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit c1ec304e2cb201466f7cab623206e7c2e7e9458b | file-specific BASE GPU-port residue with no repeatable transform
- src/Model/Solid/Finite/NeoHookean.H:36-90#3 [NUM] | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit c1ec304e2cb201466f7cab623206e7c2e7e9458b | file-specific BASE GPU-port residue with no repeatable transform
- src/Numeric/Function.H:1-6 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit 1da57132e730ebaf97014ff63bc74f0b72cbf870 | file-specific BASE GPU-port residue with no repeatable transform
- src/Numeric/Stencil.H:83-107#2 | Replaces pair return values with the device-safe Jet aggregate. | commit 757124929e04c2fcca82182bf0e3c2895867f0ad | file-specific BASE GPU-port residue with no repeatable transform
- src/Numeric/Stencil.H:785-821#3 | Replaces split-gradient pair return values with the device-safe Jet aggregate. | commit adce42a4e9e105f450bbddd4298805d1febdd346 | file-specific BASE GPU-port residue with no repeatable transform
- src/Operator/Elastic.H:108-115 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit adce42a4e9e105f450bbddd4298805d1febdd346 | file-specific BASE GPU-port residue with no repeatable transform
- src/Set/Base.H:7-12 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit 4edcfa7c095479b07fcaddb9f8160326f1cdc2d3 | file-specific BASE GPU-port residue with no repeatable transform
- src/Set/Base.H:16-21 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit bb4afc9efeb774df2af978989768235511ac6ae2 | file-specific BASE GPU-port residue with no repeatable transform
- src/Set/Matrix4_Major.H:54-93#2 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit 9470889b14f10a902dab6dd9573deefc09972e8d | file-specific BASE GPU-port residue with no repeatable transform
- src/Util/Util.H:15-20 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit dd4054056aa5ef24bde948ac1b04013b5b106eb6 | file-specific BASE GPU-port residue with no repeatable transform
- src/Util/Util.cpp:28-33 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da | file-specific BASE GPU-port residue with no repeatable transform
- src/Util/Util.cpp:138-143 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit 8ae292c331d88d40a6031b7cc2a131d1e03e98c9 | file-specific BASE GPU-port residue with no repeatable transform
- src/alamo_gpu.cc:0-0 | GPU/build compatibility support (device annotations, closure-safe types, include, or launcher plumbing). | commit c1ec304e2cb201466f7cab623206e7c2e7e9458b | file-specific BASE GPU-port residue with no repeatable transform
