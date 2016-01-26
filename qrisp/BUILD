# The main bazel configuration file for QRisp. This BUILD file defines all
# library, binary and test creation rules.

cc_library(
    name = "dataset-utils",
    srcs = ["dataset-utils.cc"],
    hdrs = ["dataset-utils.h"],
    deps = [
        "//qrisp:structure",
        "//qrisp:utils",
        "//qrisp/proto:config_proto",
        "//qrisp/proto:parameters_proto",
        "//qrisp/proto:structure_proto",
    ],
)

cc_library(
    name = "structure",
    srcs = ["structure.cc"],
    hdrs = ["structure.h"],
    deps = [
        "//qrisp:utils",
        "//qrisp/proto:structure_proto",
    ],
)

cc_library(
    name = "performance",
    srcs = ["performance.cc"],
    hdrs = ["performance.h"],
    deps = [
        "//qrisp:cluster",
        "//qrisp:dataset-utils",
        "//qrisp:structure",
        "//qrisp:utils",
        "//qrisp/proto:config_proto",
    ],
)

cc_library(
    name = "model",
    srcs = ["model.cc"],
    hdrs = ["model.h"],
    deps = [
        "//qrisp:plif",
        "//qrisp:structure",
        "//qrisp:utils",
        "//qrisp/proto:parameters_proto",
        "//third_party/googleflags:gflags",
    ],
)

cc_library(
    name = "utils",
    srcs = ["utils.h"],
    hdrs = ["utils.h"],
    deps = [
        "//third_party/googlelog:glog",
    ],
)

cc_library(
    name = "plif",
    srcs = ["plif.cc"],
    hdrs = ["plif.h"],
    deps = [
        "//qrisp:utils",
    ],
)

cc_library(
    name = "recurrences",
    srcs = [
        "recurrences.cc",
        "recurrences-nbest.cc",
    ],
    hdrs = [
        "recurrences.h",
        "recurrences-nbest.h",
    ],
    deps = [
        "//qrisp:model",
        "//qrisp:utils",
    ],
)

cc_library(
    name = "cluster",
    srcs = ["cluster.cc"],
    hdrs = ["cluster.h"],
    deps = [
        "//third_party/googlelog:glog",
        "//qrisp:utils",
    ],
)

cc_library(
    name = "qrisp-lib",
    srcs = ["qrisp-lib.cc"],
    hdrs = ["qrisp-lib.h"],
    deps = [
        "//qrisp:cluster",
        "//qrisp:dataset-utils",
        "//qrisp:learning-utils",
        "//qrisp:model",
        "//qrisp:performance",
        "//qrisp:recurrences",
        "//qrisp:sgd",
        "//qrisp:structure",
        "//qrisp:utils",
        "//qrisp/proto:parameters_proto",
        "//qrisp/proto:structure_proto",
    ],
)

# Machine learning related libaries.

cc_library(
    name = "learning-utils",
    srcs = ["learning-utils.cc"],
    hdrs = ["learning-utils.h"],
    deps = [
        "//qrisp:utils",
    ],
)

cc_library(
    name = "sgd",
    srcs = ["sgd.cc"],
    hdrs = ["sgd.h"],
    deps = [
        "//qrisp:learning-utils",
    ],
)

#
# The main binary of the project.
#

cc_binary(
    name = "qrisp",
    srcs = ["qrisp.cc"],
    linkopts = [
        "-lm",
    ],
    deps = [
        "//qrisp:dataset-utils",
        "//qrisp:qrisp-lib",
        "//qrisp:utils",
        "//third_party/googleflags:gflags",
        "//third_party/googlelog:glog",
    ],
)

#
# Various unit tests for the above defined libraries.
#

cc_test(
    name = "structure_test",
    srcs = ["structure_test.cc"],
    linkopts = [
        "-lpthread",
    ],
    deps = [
        "//qrisp:structure",
        "//qrisp:utils",
        "//third_party/googletest",
    ],
)

cc_test(
    name = "plif_test",
    srcs = ["plif_test.cc"],
    linkopts = [
        "-lpthread",
    ],
    deps = [
        "//qrisp:plif",
        "//qrisp:utils",
        "//third_party/googletest",
    ],
)

cc_test(
    name = "recurrences_test",
    srcs = ["recurrences_test.cc"],
    linkopts = [
        "-lpthread",
    ],
    deps = [
        "//qrisp:learning-utils",
        "//qrisp:recurrences",
        "//third_party/googletest",
    ],
)

cc_test(
    name = "model_test",
    srcs = ["model_test.cc"],
    linkopts = [
        "-lpthread",
        "-lm",
    ],
    deps = [
        "//qrisp:model",
        "//third_party/googletest",
    ],
)

#cc_library(
#    name = "dataset-utils",
#    srcs = ["dataset-utils.cc"],
#    hdrs = ["dataset-utils.h"],
#    deps = [
#        "//qrisp:rna-structure",
#        "//qrisp:utils",
#        "//qrisp/proto:parameters_proto",
#        "//qrisp/proto:config_proto",
#        "//qrisp/proto:structure_proto",
#    ],
#)
