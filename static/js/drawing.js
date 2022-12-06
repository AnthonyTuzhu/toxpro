
function heatmap(data) {
    // plots a heatmap of all the models in the dataset
        var margin = {top: 30, right: 30, bottom: 30, left: 30},
      width = 600 - margin.left - margin.right,
      height = 600 - margin.top - margin.bottom;

        d3.select("#stats-graph").select("svg").remove();

// append the svg object to the body of the page
    var g = d3.select("#stats-graph")
        .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
        .append("g")
      .attr("transform",
            "translate(" + margin.left + "," + margin.top + ")");

    var metrics = data[0].metrics;

    var datasets = [];

    for (var i = 0; i<data.length; i++) {
        datasets.push(data[i].name)
    }

    var x = d3.scaleBand()
      .range([ 0, width ])
      .domain(metrics)
      .padding(0.01);
        g.append("g")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(x));

        var y = d3.scaleBand()
  .range([ height, 0 ])
  .domain(myVars)
  .padding(0.01);
  g.append("g")
  .call(d3.axisLeft(y));


  var myColor = d3.scaleLinear()
  .range(["blue", "red"])
  .domain([0,1]);

  g.selectAll()
    .data(data)
    .enter()
    .append("rect")
      .attr("x", function(d, i) { return x(d) })
      .attr("y", function(d, i) { return y(d.variable) })
      .attr("width", x.bandwidth() )
      .attr("height", y.bandwidth() )
      .style("fill", function(d) { return myColor(d.value)} )


}


function plotBar(actives, inactives) {
    // bar plot that shows the active:inactive ratio
    // of a dataset

    let margin = {top: 20, right: 20, bottom: 30, left: 40};
    let w = 300 - margin.left - margin.right;
    let h = 300 - margin.top - margin.bottom;


    let svg = d3.select("#db-bar").attr("height", 300)
                .attr("width", 600);

    svg.selectAll('g').remove();

    let chart = svg.append('g').attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let yScale = d3.scaleLinear().domain([0, Math.max(actives, inactives)+10]).range([h, 0]);
    let xScale = d3.scaleBand().domain(['Actives', 'Inactives']).range([10, w]);

    let yAxis = chart.append("g")
                    .call(d3.axisLeft(yScale));
    let xAxis = chart.append("g")
                    .call(d3.axisTop(xScale));


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
}


function plotFiveFoldStats(data) {
    // plots the data from a five fold cross validation result
    // should also have data.name and data.alg_name which give the dataset name and the algorithm
    // data.results and data.metrics are arrays and
    // correspond to the actual range-scaled results and metric names, respectively
    // data should be a JSON object where the length is equal to the metrics

    let margin = {top: 20, right: 20, bottom: 30, left: 40};
    let w = 300 - margin.left - margin.right;
    let h = 300 - margin.top - margin.bottom;

    let svg = d3.select("#stats-graph").attr("height", 300)
      .attr("width", 600);

    svg.selectAll('g').remove();

    let chart = svg.append('g').attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    let yScale = d3.scaleLinear().domain([0, 1]).range([h, 0]);
    let xScale = d3.scaleBand().domain(data.metrics).range([0, w]);

    let yAxis = chart.append("g")
                    .attr("class", "yAxis")
                    .call(d3.axisLeft(yScale));
    let xAxis = chart.append("g")
                    .attr("class", "xAxis")
                    .attr("transform", "translate(0," + h + ")")
                    .call(d3.axisBottom(xScale));

    chart.selectAll("results")
      .data(data.results)
      .enter().append("rect")
      .attr("width", xScale.bandwidth())
      .attr("height", function (d) {return h-yScale(d);})
      .attr("x", function (d, i) {return xScale(data.metrics[i]);})
      .attr("y", function (d) {return yScale(d);})
      .style("fill", "black")
      .style("stroke", "black")
      .style("fill-opacity", 0.2).style("stroke-width", 3);

}