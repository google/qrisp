cc_library(
    name = "sgd",
    srcs = ["sgd.cc"],
    hdrs = ["sgd.h"],
)

cc_library(
    name = "dataset-utils",
    srcs = ["dataset-utils.cc"],
    hdrs = ["dataset-utils.h"],
)

cc_library(
    name = "learning-utils",
    srcs = ["learning-utils.cc"],
    hdrs = ["learning-utils.h"],
)

cc_library(
    name = "rna-structure",
    srcs = ["rna-structure.cc"],
    hdrs = ["rna-structure.h"],
    deps = [
        "//qrisp/proto:structure_proto",
    ],
)

cc_library(
    name = "performance",
    srcs = ["performance.cc"],
    hdrs = ["performance.h"],
)

cc_library(
    name = "model",
    srcs = ["model.cc"],
    hdrs = ["model.h"],
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
)

cc_library(
    name = "cluster",
    srcs = ["cluster.cc"],
    hdrs = ["cluster.h"],
)

#
# The main binary of the project.
#

cc_binary(
    name = "qrisp",
    srcs = ["qrisp.cc"],
    deps = [
        ":cluster",
        ":recurrences",
    ],
)
