#
#
#

cc_library(
    name = "dataset-utils",
    srcs = ["dataset-utils.cc"],
    hdrs = ["dataset-utils.h"],
    deps = [
        "//qrisp:rna-structure",
        "//qrisp:utils",
        "//qrisp/proto:parameters_proto",
        "//qrisp/proto:structure_proto",
    ],
)

cc_library(
    name = "rna-structure",
    srcs = ["rna-structure.cc"],
    hdrs = ["rna-structure.h"],
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
        "//qrisp:rna-structure",
        "//qrisp:utils",
        "//qrisp/proto:config_proto",
    ],
)

cc_library(
    name = "model",
    srcs = ["model.cc"],
    hdrs = ["model.h"],
    deps = [
        "//qrisp:rna-structure",
        "//qrisp:utils",
        "//qrisp/proto:parameters_proto",
    ],
)

cc_library(
    name = "utils",
    srcs = ["utils.cc"],
    hdrs = ["utils.h"],
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
    deps = [
        ":dataset-utils",
        ":recurrences",
        ":utils",
    ],
)
