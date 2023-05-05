
// function saveData(name, data){
//   var xhr = new XMLHttpRequest();
//   xhr.open('POST', 'write_data.php'); // 'write_data.php' is the path to the php file described above.
//   xhr.setRequestHeader('Content-Type', 'application/json');
//   xhr.send(JSON.stringify({filedata: data}));
// }

/* initialize jsPsych */
var jsPsych = initJsPsych({
  on_finish: function() {
    //jsPsych.data.displayData();
    sub_ID = jsPsych.data.get()["trials"][1]["response"]["Q0"]
    jsPsych.data.get().localSave('csv',`${sub_ID}_data.csv`);
  }
});

/* create timeline */
var timeline = [];

var count = 0;

var compounds = {
  "Air": `<img width="50" height="50" src="img/Antinomy.png">`,
  "Earth": `<img width="50" height="50" src="img/Bismuth.png">`,
  "Water": `<img width="50" height="50" src="img/Sulfur.png">`,
  "Fire": `<img width="50" height="50" src="img/Arsenic.png">`,
  "Air Earth": `<img width="50" height="50" src="img/Lead.png">`,
  "Air Earth Water": `<img width="50" height="50" src="img/Silver.png">`,
  "Earth Fire Water": `<img width="50" height="50" src="img/Gold.png">`
};

var symbols = {
  "Air":  "(\u{1F701})",
  "Earth": "(\u{1F703})",
  "Water":  "(\u{1F704})",
  "Fire": "(\u{1F702})"
}

var numCompounds = Object.keys(compounds).length;

//let tested = [["Air"], ["Earth"], ["Water"], ["Fire"]];

let tested = []

/* preload images */
var preload = {
  type: jsPsychPreload,
  images: ['img/Antinomy.png', 'img/Arsenic.png', 'img/Bismuth.png', 'img/Gold.png', 'img/Lead.png', 'img/Silver.png', 'img/Sulfur.png']
};
timeline.push(preload);

/* define welcome message trial */
var welcome = {
  type: jsPsychSurveyText,
  questions: [
    {prompt: 'Welcome! <br>The experimenter will enter your participant ID.',
    required: true,}
  ]
  //stimulus: "<p> Welcome! Please enter your participant ID</p>",
  //choices: ["Continue"]
};
timeline.push(welcome);

/* define instructions trial */
var instructionsA = {
  type: jsPsychHtmlButtonResponse,
  stimulus: `
    <p> <b> You are an alchemist. </b></p>
    <p>Before you are four basic elements: <br>
    <em>Air  </em>(\u{1F701}), <em> Earth  </em> (\u{1F703}), <em> Fire </em> (\u{1F702}) , and <em> Water </em> (\u{1F704}) </p>
    <p>In each trial of this experiment, you can choose to <br>
    add 2 to 4 of these elements into a universal solvent (\u{2192}) <br>
    to test what the compound produces. </p>

    <p>Examples:  <br>
    Air (\u{1F701}) + Water (\u{1F704}) \u{2192} ? <br> Air (\u{1F701}) + Earth (\u{1F703}) + Fire (\u{1F702}) + Water (\u{1F704}) \u{2192} ? <br>
    <div style='width: 700px;'>
  `,
  choices: ["Continue"]
};
timeline.push(instructionsA);

/* define instructions trial */
var instructionsB = {
  type: jsPsychHtmlButtonResponse,
  stimulus: `
    <center> <p> You already know the following compounds </p>
    <table> <tr> <td style= "text-align: right"> Air (\u{1F701}) \u{2192} </td> <td> <img width="50" height="50" src="img/Antinomy.png"> </td> </tr>
    <tr> <td style= "text-align: right"> Earth (\u{1F703}) \u{2192} </td> <td> <img width="50" height="50" src="img/Bismuth.png"> </td> </tr>
    <tr> <td style= "text-align: right"> Fire (\u{1F702}) \u{2192} </td> <td> <img width="50" height="50" src="img/Arsenic.png"> </td> </tr>
    <tr> <td style= "text-align: right"> Water (\u{1F704}) \u{2192} </td> <td> <img width="50" height="50" src="img/Sulfur.png"> </td> </tr> </table> </center>
    
  `,
  choices: ["Continue"]
};
timeline.push(instructionsB);


var instructionsC = {
  type: jsPsychHtmlButtonResponse,
  stimulus: `
    <p> Given this knowledge, mathematically, there are 11 possible ways to combine each of the four elements with the solvent. <br>
    There is a one-to-one relationship between compounds and their products <br>
    (so there's no need to test the same compound more than once). <br>
    Order also does not matter. <br>
    <p> It's your job to figure out what each compound produces. <br>
    Due to a shortage of solvent, <em> you can only test 8 compounds, </em>but <b> you should try to finish in as few trials as possible. </b> <br>
    After 8 trials, you will be required to report what you think each compound produces. </p>
    <p> If you think you know all of the compounds before your 8 trials are completed, click the "I'm done" button. </p>
    <p> When you're ready to start testing, click continue. </p>
    <div style='width: 700px;'>
  `,
  choices: ["Continue"]
};
timeline.push(instructionsC);


// Split this up more
// 10 might be too few
// You should do as few as possible
/* define instructions trial */

var trial = {
  type: jsPsychSurveyMultiSelect,
  questions: [
    {
      prompt: "<center> <p> What compound would you like to test? <b>(Select 2, 3, or 4)</b> </p> </center>", 
      options: ["Air (\u{1F701})", "Earth (\u{1F703})", "Fire (\u{1F702})", "Water (\u{1F704})"], 
      horizontal: true,
      required: true,
      name: 'Trial'
    }, 
    {
      prompt: function(allTested=tested, emergentCompounds = compounds, N = numCompounds, symbol_list = symbols) {
        known_compounds = `<center> <table style="width:100%"> 
          <tr> 
            <td> <center> Air (\u{1F701}) </center></td>
            <td> <center> Earth (\u{1F703}) </center></td>
            <td> <center> Fire (\u{1F702}) </center></td>
            <td> <center> Water (\u{1F704}) </center></td>
          </tr>
          <tr> 
            <td> <center> \u{2193} </center></td>
            <td> <center> \u{2193} </center></td>
            <td> <center> \u{2193} </center></td>
            <td> <center> \u{2193} </center></td>
          </tr>
          <tr> 
            <td> <center> <img width="50" height="50" src="img/Antinomy.png"> </center></td>
            <td> <center> <img width="50" height="50" src="img/Bismuth.png"> </center></td>
            <td> <center> <img width="50" height="50" src="img/Arsenic.png"> </center></td>
            <td> <center> <img width="50" height="50" src="img/Sulfur.png"> </center></td>
          </tr> </table> </center>`

        if (allTested.length != 0){
          var all_phrases = []

          for (let j = 0; j < allTested.length; j++){
            var output = []
            var words = []
            for (let i = 0; i < N; i++) {

              for (let k = 0; k < allTested[j].length; k++){
                words.push(allTested[j][k].split(" ")[0])
              }

              key_as_list = Object.keys(emergentCompounds)[i].split(" ")

              //console.log(allTested[j])
              if (key_as_list.every(val => words.includes(val))){
                output.push(emergentCompounds[Object.keys(emergentCompounds)[i]])
              }
            }
          
            input_string = allTested[j].join(' + ')
            output_string = output.join(' ')



            all_phrases.push(`<tr> <td style= "text-align: right"> ${j+1}: ${input_string} \u{2192} </td> <td> ${output_string} </td> </tr>`)


          }

          
          return `<hr> <center> Known Compounds: </center> ${known_compounds} <hr> <p> <center> Tested Compounds: </p> <table> ${all_phrases.join(' ')} </table> </center>`
        }
        else{
          return `<hr> <center> Known Compounds: </center> ${known_compounds} <hr> <center>  Tested Compounds: </p> No compounds have been tested yet </center>`
        }
      },
      options:[]
    }
  ]
};

//timeline.push(trial);


var output_trial = {
  type: jsPsychHtmlButtonResponse,
  stimulus: function(allTested=tested, emergentCompounds = compounds, N = numCompounds) {
    var words = []
    var LastTrialResponse = jsPsych.data.getLastTrialData()["trials"][0]["response"]["Trial"]//.filter({trials: 0});
    console.log(LastTrialResponse)

    for (let k = 0; k < LastTrialResponse.length; k++){
      words.push(LastTrialResponse[k].split(" ")[0])
    }
    //


    allTested.push(LastTrialResponse)

    var output = []

    for (let i = 0; i < N; i++) {
      key_as_list = Object.keys(emergentCompounds)[i].split(" ")

      if (key_as_list.every(val => words.includes(val))){
        output.push(emergentCompounds[Object.keys(emergentCompounds)[i]])
      }
    }

    input_string = words.join(' + ')
    output_string = output.join(' ')


    return `<table> <tr> <td> ${input_string} \u{2192} </td> <td> ${output_string} </td> </tr> </table>`;

  },
  choices: ["Continue", "I'm done"]

}


var loop_node = {
    timeline: [trial, output_trial],
    loop_function: function(data, c = count){
      //console.log(c)
        if(c > 6 || data["trials"][1]["response"] == 1){
            return false;
        } else {
          count = count+1
            return true;
        }
    }
}

timeline.push(loop_node);

var end = {
  type: jsPsychHtmlButtonResponse,
  stimulus: function(allTested=tested, emergentCompounds = compounds, N = numCompounds, symbol_list = symbols) {
        var all_phrases = []
        known_compounds = `<center> <table style="width:75%"> 
          <tr> 
            <td> <center> Air (\u{1F701}) </center></td>
            <td> <center> Earth (\u{1F703}) </center></td>
            <td> <center> Fire (\u{1F702}) </center></td>
            <td> <center> Water (\u{1F704}) </center></td>
          </tr>
          <tr> 
            <td> <center> \u{2193} </center></td>
            <td> <center> \u{2193} </center></td>
            <td> <center> \u{2193} </center></td>
            <td> <center> \u{2193} </center></td>
          </tr>
          <tr> 
            <td> <center> <img width="50" height="50" src="img/Antinomy.png"> </center></td>
            <td> <center> <img width="50" height="50" src="img/Bismuth.png"> </center></td>
            <td> <center> <img width="50" height="50" src="img/Arsenic.png"> </center></td>
            <td> <center> <img width="50" height="50" src="img/Sulfur.png"> </center></td>
          </tr> </table> </center>`

        for (let j = 0; j < allTested.length; j++){
          var output = []
          var words = []
          for (let i = 0; i < N; i++) {

            for (let k = 0; k < allTested[j].length; k++){
              words.push(allTested[j][k].split(" ")[0])
            }

            key_as_list = Object.keys(emergentCompounds)[i].split(" ")

            //console.log(allTested[j])
            if (key_as_list.every(val => words.includes(val))){
              output.push(emergentCompounds[Object.keys(emergentCompounds)[i]])
            }
          }
        
          input_string = allTested[j].join(' + ')
          output_string = output.join(' ')



          all_phrases.push(`<tr> <td style= "text-align: right"> ${j+1}: ${input_string} \u{2192} </td> <td> ${output_string} </td> </tr>`)


        }
        sub_ID = jsPsych.data.get()["trials"][1]["response"]["Q0"]
        jsPsych.data.get().localSave('csv',`${sub_ID}_data.csv`)

        return `<p> Congratulations! <br> The experimenter will now ask you to record what you think each compound produces. <br>
        For reference, below are the compounds you tested. </p>
        <center> <hr> Known Compounds: <br> ${known_compounds} <hr>
        <center> Tested Compounds: </p> <table> ${all_phrases.join(' ')} </table> </center>`
      },

  //"<p> Congratulations! The experimenter will instruct you on what to do next</p>",
  choices: []
};
timeline.push(end);



/* start the experiment */
//jsPsych.run([timeline, loop_node]);
jsPsych.run(timeline);
