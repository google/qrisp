package(default_visibility = ["//visibility:public"])

load("/tools/build_rules/proto", "proto_library")

licenses(["notice"])  # Apache2

filegroup(
    name = "public",
    srcs = [
        "config.proto",
        "parameters.proto",
        "structure.proto",
    ],
)

# Proto storing relevant stats from logs for a given criteria.
proto_library(
    name = "structure_proto",
    src = "structure.proto",
)

# Proto storing all settings relevant for learning.
proto_library(
    name = "config_proto",
    src = "config.proto",
)

# Proto storing relevant stats from logs for a given criteria.
proto_library(
    name = "parameters_proto",
    src = "parameters.proto",
    deps = [":config_proto"],
)
