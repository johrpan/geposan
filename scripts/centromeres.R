library(data.table)

# This depends on the current development version of this package to be
# installed in order to retrieve the distance data.
distances <- geposan::distances

setorder(distances, "start_position")

chromosomes <- distances[,
    .(chromosome_length = max(end_position)),
    by = c("species", "chromosome_name")
]

handle_chromosome <- function(species_, chromosome_name_) {
    chromosome_distances <- distances[
        species == species_ & chromosome_name == chromosome_name_
    ]

    chromosome_distances[,
        center_position := mean(c(start_position, end_position)),
        by = gene
    ]

    setorder(chromosome_distances, center_position)

    chromosome_distances[1, gap_size := 0]

    for (index in 2:nrow(chromosome_distances)) {
        previous_end_position <- chromosome_distances[index - 1, end_position]
        chromosome_distances[
            index,
            gap_size := start_position - previous_end_position
        ]
    }

    chromosome_gaps <- chromosome_distances[order(gap_size, decreasing = TRUE)]

    # Retrieve the two largest gap sizes.
    max_gap_size <- chromosome_gaps[1, gap_size]
    next_gap_size <- chromosome_gaps[2, gap_size]

    list(
        centromere_start = chromosome_distances[
            which.max(gap_size) - 1,
            end_position
        ],
        centromere_end = chromosome_distances[
            which.max(gap_size),
            start_position
        ],
        centromere_dominance = max_gap_size / next_gap_size
    )
}

chromosomes[,
    c("centromere_start", "centromere_end", "centromere_dominance") :=
        handle_chromosome(species, chromosome_name),
    by = c("species", "chromosome_name")
]
