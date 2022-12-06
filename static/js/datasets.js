function plotBar(actives, inactives) {
    // bar plot that shows the active:inactive ratio
    // of a dataset



    var cardBody = d3.select("#dataset-bar-body");
    currentWidth = parseInt(cardBody.style("width"));

    var margin = {top: 20, right: 20, bottom: 30, left: 40};
    var w = currentWidth - margin.left - margin.right;
    var h = 400 - margin.top - margin.bottom;

    var svg = d3.select("#dataset-bar")
        .attr("width", w + margin.left + margin.right)
        .attr("height", h + margin.top + margin.bottom);

    svg.selectAll('g').remove();

    var chart = svg.append('g')
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var xScale = d3.scaleBand()
        .domain(['Actives', 'Inactives'])
        .range([margin.left, w - margin.right])
        .paddingInner(0.10)
        .paddingOuter(0.30);


    var yScale = d3.scaleLinear().domain([0, Math.max(actives, inactives)+10]).range([h, 0]);





    chart.selectAll("act")
      .data([actives])
      .enter().append("rect")
      .attr("width", xScale.bandwidth())
      .attr("height", function (d) {
          return h-yScale(d);
      })
      .attr("x", xScale('Actives'))
      .attr("y", function (d) {
          return yScale(d);
      })
      .style("fill", "red")
      .style("stroke", "black")
      .style("fill-opacity", 0.4).style("stroke-width", 3);


    chart.selectAll("inact")
      .data([inactives])
      .enter().append("rect")
      .attr("width", xScale.bandwidth())
      .attr("height", function (d) {
          return h-yScale(d);
      })
      .attr("x", xScale('Inactives'))
      .attr("y", function (d) {
          return yScale(d)
      })
      .style("fill", "green")
      .style("stroke", "black")
      .style("fill-opacity", 0.4).style("stroke-width", 3);


    var yAxis = chart.append("g")
                    .attr("transform", "translate(" + margin.left + ",0)")
                    .call(d3.axisLeft(yScale))
                    .selectAll("text")
                    .style("font-size", 18);

    var xAxis = chart.append("g")
                    .attr("transform", "translate(0," + yScale(0) + ")")
                    .call(d3.axisBottom(xScale))
                    .selectAll("text")
                    .style("font-size", 30);


}

function updateOverview(data) {

    $("#this_ds_text").text(data.name);
    $("#num_act_text").text(data.actives);
    $("#num_inact_text").text(data.inactives);
    $("#tot_text").text(data.set_type);

}

function updateCompoundTable(data) {

    var tableScroll = $("#compound-overview").empty();
    var table = $('<table></table>').addClass("table table-bordered table-striped mb-0");
    tableScroll.append(table);

    // table header

        var head = $("                            <thead>\n" +
        "                              <tr>\n" +
        "                                  <th scope=\"row\" >Identifier</th>\n" +
        "                                  <th scope=\"row\" >Activity</th>\n" +
        "                              </tr>\n" +
        "                            </thead>");

    table.append(head);

    var tableBody = $('<tbody></tbody>')
    table.append(tableBody);

    for (var i = 0; i < data.compounds.length; i++) {


        var row = $("<tr></tr>");

        tableBody.append(row);

        var rowHead = $("<th scope=\"row\"></th>");

        rowHead.text(data.compounds[i]);
        row.append(rowHead);

        var activity = $("<td></td>").text(data.activities[i].toString());

        // color cell according to activity
        if (data.activities[i] == 1) {
            activity.css("background-color", "rgba(255, 0, 0, 0.3)");
        } else {
            activity.css("background-color", "rgba(0, 255, 0, 0.3)");
        }
        row.append(activity);

    }

    // set height to be same as the activity overview

    tableScroll.attr("height")

}

function updateDataset() {
    var e = document.getElementById("dataset-selection");
    var currentProfile = e.options[e.selectedIndex].value;

    var queryUrl = $SCRIPT_ROOT + "get_dataset_overview/" + currentProfile;
    var dataset_data = JSON.parse(getResponseFromURL(queryUrl));

    updateOverview(dataset_data);
    plotBar(dataset_data.actives, dataset_data.inactives);
    updateCompoundTable(dataset_data);

}

function deleteDataset() {
    var e = document.getElementById("dataset-selection");
    var currentDataset = e.options[e.selectedIndex].value;


    data = {
        dataset_name: currentDataset
    }

    postData('/delete_dataset', data);

    location.reload();
}

function addInhouseDataset(dsName) {


    data = {
        inhouse_ds_name: dsName
    }
    console.log(data);
    postData('/add_inhouse_dataset', data);

    location.reload();
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
