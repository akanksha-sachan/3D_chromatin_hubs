#!/usr/bin/env node
import Straw from "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/hic-straw/src/straw.js"
import NodeLocalFile from "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/hic-straw/src/io/nodeLocalFile.mjs"
const path = "/Users/Akanksha/MaGroup/Genomic Hubs/workflow/data/GM12878/ENCODE_hg38/ENCFF318GOM.hic"
const nodeLocalFile = new NodeLocalFile({path})
const straw = new Straw({file: nodeLocalFile})

// Get the normalization options as an array
straw.getNormalizationOptions()
    .then(function (normOptions) {
        console.log(JSON.stringify({ normalizationOptions: normOptions }));
    })
    .catch(error => {
        console.error(error);
        process.exit(1);
    });