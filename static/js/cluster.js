function getRandomColor() {
  var letters = '0123456789ABCDEF';
  var color = '#';
  for (var i = 0; i < 6; i++) {
    color += letters[Math.floor(Math.random() * 16)];
  }
  return color;
}



function changeClusteringName() {
    selector = $("#profile-selection");
    currentProfile = selector.find(":selected").text().trim();
    clusteringFilename = $("#clustering_filename");
    clusteringFilename.val(currentProfile + "_clustering");

};

function updateMaxNoClusters(maxNumber) {
   $("#n_clusters").attr("max", maxNumber);
};


function refreshProfileClusterPage() {

    // replaces all the profile data on the cluster page

    var currentProfile = $("#profile-selection").find(":selected").text().trim();

    var queryUrl = $SCRIPT_ROOT + "get_bioprofile/" + currentProfile;
    var profile_data = JSON.parse(getResponseFromURL(queryUrl));

    var tsElement = document.getElementById('this_ts');
    tsElement.innerHTML = profile_data.meta.training_set;

    var tsElement = document.getElementById('cmps');
    tsElement.innerHTML = profile_data.meta.num_cmps;

    var tsElement = document.getElementById('aids');
    tsElement.innerHTML = profile_data.meta.num_aids;

    var tsElement = document.getElementById('tot_acts');
    tsElement.innerHTML = profile_data.meta.num_total_actives;

    var tsElement = document.getElementById('tot_inacts');
    tsElement.innerHTML = profile_data.meta.num_total_inactives;

    updateMaxNoClusters(profile_data.meta.num_aids);
}






function networkGraph (graph) {


    d3.select("#cluster-graph").select("svg").remove();

    var cardBody = d3.select("#cluster-graph-card-body");

    var body = d3.select("#cluster-graph"),
        width = parseInt(cardBody.style("width"))-100,
        height = parseInt(cardBody.style("width"))-100,
        radius = 6;

    // cluster colors
    // fill = d3.scaleOrdinal(d3.schemeCategory20).domain([0, graph.nodes.length-1]);

    fill = [];
    for (i = 0; i < graph.nodes.length; i++) {
            fill.push(getRandomColor());
        }


    var svg = body.append("svg");

    svg.attr("width", width);
    svg.attr("height", height);

    var simulation = d3.forceSimulation()
        .force("link", d3.forceLink().id(function(d) { return d.id; })
                                     .distance(100)
                                     // .strength(function(d) { return 1-d.weight; })
        )
        .force("charge", d3.forceManyBody())
        .force("center", d3.forceCenter(width / 2, height / 2));



      var link = svg.append("g")
                    .style("stroke", "#aaa")
                    .selectAll("line")
                    .data(graph.links)
                    .enter().append("line");

      var node = svg.append("g")
                .attr("class", "nodes")
      .selectAll("circle")
                .data(graph.nodes)
      .enter().append("circle")
              .attr("r", radius)
              .style("fill", function(d) { console.log(d.class);
                                          console.log(fill[d.class]);
                                            return fill[d.class]; })
             .style("stroke", "#969696")
             .style("stroke-width", "1px")
             .attr("classLabel", function(d) { return d.class; })
             .attr("id", function(d) { return "AID" + d.id; })
             .attr("AID", function(d) { return d.id; })
              .call(d3.drag()
                  .on("start", dragstarted)
                  .on("drag", dragged)
                  .on("end", dragended)
              );

      var label = svg.append("g")
          .attr("class", "labels")
          .selectAll("text")
          .data(graph.nodes)
          .enter().append("text")
            .attr("class", "label")
            .text(function(d) { return d.name; });

      simulation
          .nodes(graph.nodes)
          .on("tick", ticked);

      simulation.force("link")
          .links(graph.links);

      function ticked() {
        link
            .attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });

        node
             .attr("cx", function (d) { return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
             .attr("cy",  function(d) { return d.y = Math.max(radius, Math.min(height - radius, d.y)); });

        label
                .attr("x", function(d) { return d.x; })
                .attr("y", function (d) { return d.y; })
                .style("font-size", "20px").style("fill", "rgba(0, 0, 0, 0.4)");
      }


        function dragstarted(d) {
          if (!d3.event.active) simulation.alphaTarget(0.3).restart();
          d.fx = d.x;
          d.fy = d.y;
        }

        function dragged(d) {
          d.fx = d3.event.x;
          d.fy = d3.event.y;
        }

        function dragended(d) {
          if (!d3.event.active) simulation.alphaTarget(0);
          d.fx = null;
          d.fy = null;
        }
}



function findMinIdx(Z, n, length) {
    if (Z[n][0] >= length) {
        return Math.min(findMinIdx(Z, Z[n][0]-length, length))
    } else if (Z[n][1] >= length) {
        return Math.min(Z[n][0], findMinIdx(Z, Z[n][1]-length, length))
    } else {
        return Math.min(Z[n][0], Z[n][1])
    }
}


function findAssocIdx(Z, n, length) {
    if ((Z[n][1] >= length) && (Z[n][0] >= length)) {
        return findAssocIdx(Z, Z[n][0]-length, length).concat(findAssocIdx(Z, Z[n][1]-length, length))
    } else if (Z[n][1] >= length) {
        return [Z[n][0]].concat(findAssocIdx(Z, Z[n][1]-length, length))
    } else {
        return [Z[n][0], Z[n][1]]
    }
}

function merge(graph, n) {




    if (n >= 0) {

        var assocIdx = findAssocIdx(graph.linkage, n, graph.nodes.length);

        var prevClasses = [];
        for (index = 0; index < assocIdx.length; index++) {
            prevClasses.push(graph.nodes[assocIdx[index]].class);
        }

        var newClass = Math.min(...prevClasses);
        for (index = 0; index < assocIdx.length; index++) {
            graph.nodes[assocIdx[index]].class = newClass;
        }


    }



}


function unmerge(graph, n) {
    var cluserOneIdx = graph.linkage[n][0];
    var cluserTwoIdx = graph.linkage[n][1];
    console.log(cluserOneIdx, cluserTwoIdx)
    if (cluserOneIdx < graph.nodes.length) {
            graph.nodes[cluserOneIdx].class = graph.nodes[cluserOneIdx].index
    } else {
        cluserOneIdx = findAssocIdx(graph.linkage, cluserOneIdx-graph.nodes.length, graph.nodes.length);
        console.log(cluserOneIdx)
        // var clusterClasses = [];
        // for (i = 0; i < cluserOneIdx.length; i++) {
        //     clusterClasses.push(graph.nodes[cluserOneIdx[i]].class)
        // }

        var newClass = Math.min(...cluserOneIdx);
        console.log(newClass);
        for (i = 0; i < cluserOneIdx.length; i++) {
            graph.nodes[cluserOneIdx[i]].class = newClass;
        }
    }

    if (cluserTwoIdx < graph.nodes.length) {
            graph.nodes[cluserTwoIdx].class = graph.nodes[cluserTwoIdx].index
    } else {
        cluserTwoIdx = findAssocIdx(graph.linkage, cluserTwoIdx-graph.nodes.length, graph.nodes.length);
        console.log(cluserTwoIdx)
        // var clusterClasses = [];
        // for (i = 0; i < cluserTwoIdx.length; i++) {
        //     clusterClasses.push(graph.nodes[cluserTwoIdx[i]].class)
        // }

        var newClass = Math.min(...cluserTwoIdx);
        for (i = 0; i < cluserTwoIdx.length; i++) {
            graph.nodes[cluserTwoIdx[i]].class = newClass;
        }
    }

}


function updateColors() {



    // var fill = d3.scaleOrdinal(d3.schemeCategory20).domain([0, graph.nodes.length]);
    nodes = d3.selectAll("circle");

    nodes.data(graph.nodes);
    nodes.filter(function(d) {
      return d.isChanged;
    }).transition().duration(1000)
        .style("r", 16)
        .transition().duration(750).style("fill", function (d) {return fill[d.class]})
    .transition().duration(1000)
        .style("r", 6)




}



function postData(url, data) {
    // function that uses fetch model to send the current filters to the flask function


    fetch(url, {

            method: 'POST',

            headers: {
                'Content-type': 'application/json'
            },

            body: JSON.stringify(data)
        }
        )
}

function sendClusterData() {
    // all the circles on the HTML should be PubChem AIDs

    var results = {};
    results.clusterAssignments = [];

    for (i = 0; i < graph.nodes.length; i++) {
        var data = {};

        data.aid = +graph.nodes[i].id;
        data.classLabel = +graph.nodes[i].class;

        results.clusterAssignments.push(data);
    }

    results.currentClustering = $("#cluster-selection").find(":selected").text().trim();

    data = {
        results: results
    };
    console.log(data);
    postData('/sendClusterData', data);

}