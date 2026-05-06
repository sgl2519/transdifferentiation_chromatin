#!/usr/bin/env gawk

BEGIN {
    FS = OFS = "\t"
}

# Extract attribute value (gene_id, transcript_id)
function get_attr(attr, key, a, i, v) {
    split(attr, a, ";")
    for (i in a) {
        if (a[i] ~ key"[[:space:]]+") {
            v = a[i]
            gsub(".*" key "[[:space:]]+\"", "", v)
            gsub("\".*", "", v)
            return v
        }
    }
    return "NA"
}

# Read exon lines
$3 == "exon" {
    g = get_attr($9, "gene_id")
    t = get_attr($9, "transcript_id")

    if (g == "NA")
        next

    n[g]++

    chr[g, n[g]]    = $1
    src[g, n[g]]    = $2
    start[g, n[g]]  = $4
    end[g, n[g]]    = $5
    score[g, n[g]]  = $6
    strand[g, n[g]] = $7
    frame[g, n[g]]  = $8
    attr[g, n[g]]   = $9
    gene[g, n[g]]   = g
    tx[g, n[g]]     = t
}

END {
    for (g in n) {

        # Single-exon gene
        if (n[g] == 1) {
            print chr[g,1], src[g,1], "exon",
                  start[g,1], end[g,1],
                  score[g,1], strand[g,1], frame[g,1],
                  attr[g,1] "; exon_class \"unique\";",
                  gene[g,1], tx[g,1]
            continue
        }

        # Sort exons by strand-aware order
        for (i = 1; i <= n[g]; i++) {
            for (j = i + 1; j <= n[g]; j++) {

                if ((strand[g,1] == "+" && start[g,i] > start[g,j]) ||
                    (strand[g,1] == "-" && end[g,i]   < end[g,j])) {
                    ttmp = start[g,i];  start[g,i]  = start[g,j];  start[g,j]  = ttmp
                    ttmp = end[g,i];    end[g,i]    = end[g,j];    end[g,j]    = ttmp
                    ttmp = chr[g,i];    chr[g,i]    = chr[g,j];    chr[g,j]    = ttmp
                    ttmp = src[g,i];    src[g,i]    = src[g,j];    src[g,j]    = ttmp
                    ttmp = score[g,i];  score[g,i]  = score[g,j];  score[g,j]  = ttmp
                    ttmp = strand[g,i]; strand[g,i] = strand[g,j]; strand[g,j] = ttmp
                    ttmp = frame[g,i];  frame[g,i]  = frame[g,j];  frame[g,j]  = ttmp
                    ttmp = attr[g,i];   attr[g,i]   = attr[g,j];   attr[g,j]   = ttmp
                    ttmp = gene[g,i];   gene[g,i]   = gene[g,j];   gene[g,j]   = ttmp
                    ttmp = tx[g,i];     tx[g,i]     = tx[g,j];     tx[g,j]     = ttmp
                }
            }
        }

        # Emit output
        for (i = 1; i <= n[g]; i++) {

            if (i == 1)
                cls = "first"
            else if (i == n[g])
                cls = "last"
            else
                cls = "other"

            print chr[g,i], src[g,i], "exon",
                  start[g,i], end[g,i],
                  score[g,i], strand[g,i], frame[g,i],
                  attr[g,i] "; exon_class \"" cls "\";",
                  gene[g,i], tx[g,i]
        }
    }
}

