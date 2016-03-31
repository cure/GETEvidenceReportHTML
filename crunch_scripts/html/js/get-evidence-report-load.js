function variant_link(variant_name) {
  var base_url = "http://evidence.pgp-hms.org";
  var url = base_url + "/" + encodeURIComponent(variant_name);
  return "<a href=\"" + url + "\">" + variant_name + "</a>";
}

function populate_tables(variant_list) {
  var suff_count = 0;
  var insuff_count = 0;

  var insuff_table = $("table-insuff");

  suff_rare = [];
  suff_variants = [];
  insuff_variants = [];

  for (var i=0; i<variant_list.length; i++) {

    var variant = variant_list[i];

    if (variant.suff_eval) {

      var rare_flag = false;

      suff_count++;

      var ele = {};

      var gene = variant.gene;
      var aa_change = variant.amino_acid_change;

      // Use dbSNP if gene or aa_change is undefined
      //
      if ((typeof(gene)==="undeifned") || (typeof(aa_change)=="undefined")) {
        var dbsnps = variant.dbSNP.split(",");
        for (var ii=0; ii<dbsnps.length; ii++) {
          ele["variant_name"] = ((ii==0) ? "" : ",");
          ele["variant_name"] += variant_link(dbsnps[ii]);
        }
      } else {
        ele["variant_name"] = variant_link(variant.gene + "-" + variant.amino_acid_change);
      }

      ele["clinical_importance"] = variant.qualified_impact;
      ele["impact"] = variant.impact;

      // Determine if the alternate is homozygous or heterozygous
      //
      alts = variant.genotype.split("/");
      if (alts.length==1) {
        ele["zygosity"] = "Hom"
      } else {
        var n = 0;
        if (alts[0] != variant.ref_allele) { n++; }
        if (alts[1] != variant.ref_allele) { n++; }

        if (n==2) { ele["zygosity"] = "Hom"; }
        else { ele["zygosity"] = "Het"; }
      }
      ele["ref_allele"] = variant.ref_allele;
      ele["alt_allele"] = variant.genotype;

      // Allele frequency
      //
      if (("num" in variant) && ("denom" in variant)) {
        ele["allele_freq"] = variant.num.toString() + "/" + variant.denom.toString();

        var af = 1;
        var af_num = parseFloat(variant.num);
        var af_den = parseFloat(variant.denom);

        if (isNaN(af_num) || isNaN(af_den)) { af = 1; }
        else if (af_den>0) {
          af = af_num / af_den;

          ele["allele_freq"] += " (" + (100*af).toFixed(2).toString() + "%)";

          // If it's pathogenic and it either has allele frequency < 0.25 or is likely or well established,
          // flag it.
          //
          if ( ((af<=0.025) || /likely|well-established/i.exec(ele.clinical_importance)) &&
               (/pathogenic/i.exec(ele.clinical_importance)) ) {
            rare_flag = true;
          }
        }

      } else {
        ele["allele_freq"] = "Unknown";
      }
      ele["summary"] = variant.summary_short;

      suff_variants.push(ele);
      if (rare_flag) {
        suff_rare.push(ele);
      }


    } else {

      // There are too many entries if we consider autoscore of 0.
      //
      if (("autoscore" in variant) && (parseInt(variant.autoscore) == 0)) {
        continue;
      }

      var ele = {};
      ele["variant_name"] = variant_link(variant.gene + "-" + variant.amino_acid_change);
      ele["prioritization"] = variant.autoscore;

      // Allele frequency
      //
      if (("num" in variant) && ("denom" in variant)) {
        ele["allele_freq"] = variant.num.toString() + "/" + variant.denom.toString();

        var af = 1;
        var af_num = parseFloat(variant.num);
        var af_den = parseFloat(variant.denom);

        if (isNaN(af_num) || isNaN(af_den)) { af = 1; }
        else if (af_den>0) {
          af = af_num / af_den;

          ele["allele_freq"] += " (" + (100*af).toFixed(2).toString() + "%)";
        }

      } else {
        ele["allele_freq"] = "Unknown";
      }

      ele["num_articles"] = variant.n_articles;

      alts = variant.genotype.split("/");
      if (alts.length==1) {
        ele["zygosity"] = "Hom"
      } else {
        var n = 0;
        if (alts[0] != variant.ref_allele) { n++; }
        if (alts[1] != variant.ref_allele) { n++; }

        if (n==2) { ele["zygosity"] = "Hom"; }
        else { ele["zygosity"] = "Het"; }
      }
      ele["reason"] = variant.impact;

      ele["summary"] = "reason:" + variant.summary_short;


      if (variant_list[i].impact == "pathogenic") {
      }

      insuff_variants.push(ele);

      insuff_count++;
    }

  }

  insuff_variants.sort(function(x,y) {
    var a = parseInt(x.prioritization);
    var b = parseInt(y.prioritization);
    return b-a;
  });

  var $report_table = $("#table-report");
  $(function() {
    $report_table.bootstrapTable({data: suff_rare})
  });

  var $report_table_all = $("#table-report-all");
  $(function() {
    $report_table_all.bootstrapTable({data: suff_variants})
  });

  var $insuff_table = $("#table-insuff");
  $(function() {
    $insuff_table.bootstrapTable({data: insuff_variants});
  });

}

function allele_freq_sorter(a,b) {

  val_a = 1.0;
  val_b = 1.0;

  if ((a == "Unknown") || (a == "unknown") || (a == "Unk") || (a == "unk")) {

    // Force 'unknown' allele frequenc to be above maximum (1.0)
    //
    val_a = 1.125;

  } else {
    var a_ratio = a.split(" ")[0].split("/");
    var a_num = parseFloat(a_ratio[0]);
    var a_den = parseFloat(a_ratio[1]);
    if (a_den > 0) { val_a = a_num / a_den; }

    // clamp
    //
    if (val_a>1.0) { val_a = 1.0; }
    if (val_a<0.0) { val_a = 0.0; }

  }
  if ((b == "Unknown") || (b == "unknown") || (b == "Unk") || (b == "unk")) {

    // Force 'unknown' allele frequenc to be above maximum (1.0)
    //
    val_b = 1.125;

  } else {
    var b_ratio = b.split(" ")[0].split("/");
    var b_num = parseFloat(b_ratio[0]);
    var b_den = parseFloat(b_ratio[1]);
    if (b_den > 0) { val_b = b_num / b_den; }

    // clamp
    //
    if (val_b>1.0) { val_b = 1.0; }
    if (val_b<0.0) { val_b = 0.0; }
  }

  if (val_a < val_b) { return -1; }
  if (val_a > val_b) { return  1; }
  return 0;
}

function remove_spinner() {

  $("#spinner").css("display", "none");
  $("#main").css("display", "inline");

}

$(document).ready( function() {
  data_url = "get-evidence.json";

  // The report isn't proper JSON.  Each line
  // in the report is a JSON line, so needs to
  // be parsed on a line by line basis.  Do
  // an local AJAX request for the report
  // as a text file.
  //

  $.ajax({
    url: data_url,
    method: "GET",
    dataType: "text",
    success: function(data) {

      var variant_list = [];
      var lines = data.split("\n");
      for (var i=0; i<lines.length; i++) {
        var line = lines[i].trim();
        if (line.length==0) { continue; }
        var ele = JSON.parse(line);
        variant_list.push(ele);
      }

      populate_tables(variant_list);

      remove_spinner();

    },
    error: function(e) {
      console.log(">>> error", e);
    }
  });

});
