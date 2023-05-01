
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
    sub_ID = jsPsych.data.get()["trials"][0]["response"]["Q0"]
    jsPsych.data.get().localSave('csv',`${sub_ID}_data.csv`);
  }
});

/* create timeline */
var timeline = [];

var count = 0;

var compounds = {
  "Air": "Antimony (\u{2641})",
  "Earth": "Bismuth (\u{2646})",
  "Water": "Arsenic (\u{1F73A})",
  "Fire": "Sulfer (\u{1F70D})",
  "Air Earth": "Lead (\u{2644})",
  "Air Earth Water": "Silver (\u{263D})",
  "Earth Fire Water": "Gold (\u{2609})"
};

var symbols = {
  "Air":  "(\u{1F701})",
  "Earth": "(\u{1F703})",
  "Water":  "(\u{1F704})",
  "Fire": "(\u{1F702})"
}

var numCompounds = Object.keys(compounds).length;

let tested = [];

// /* preload images */
// var preload = {
//   type: jsPsychPreload,
//   images: ['img/blue.png', 'img/orange.png']
// };
// timeline.push(preload);

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
    add 1 to 4 of these elements into a universal solvent (\u{2192}) <br>
    to test what the compound produces. </p>

    <p>Examples: <br> Water (\u{1F704}) \u{2192} ? <br> Air (\u{1F701}) + Water (\u{1F704}) \u{2192} ? <br> Air (\u{1F701}) + Earth (\u{1F703}) + Fire (\u{1F702}) + Water (\u{1F704}) \u{2192} ? </p>
    <div style='width: 700px;'>
  `,
  choices: ["Continue"]
};
timeline.push(instructionsA);

/* define instructions trial */
var instructionsD = {
  type: jsPsychHtmlButtonResponse,
  stimulus: `
    <p> There are 16 ways to combine each of the four elements with the solvent. <br>
    There is a one-to-one relationship between compounds and their products <br>
    (so there's no need to test the same compound more than once). <br>
    <p> It's your job to figure out what each compound produces. <br>
    Due to a shortage of solvent, <em> you can only test 10 compounds. </em><br>
    After 10 trials, you will be required to report what you think each compound produces. </p>
    <p> If you think you know all of the compounds before your 10 trials are completed, click the "I'm done" button. </p>
    <p> When you're ready to start testing, click continue. </p>
    <div style='width: 700px;'>
  `,
  choices: ["Continue"]
};
timeline.push(instructionsD);

var trial = {
  type: jsPsychSurveyMultiSelect,
  questions: [
    {
      prompt: "What compound would you like to test?", 
      options: ["Air (\u{1F701})", "Earth (\u{1F703})", "Water (\u{1F704})", "Fire (\u{1F702})"], 
      horizontal: true,
      required: true,
      name: 'Trial'
    }, 
    {
      prompt: function(allTested=tested, emergentCompounds = compounds, N = numCompounds, symbol_list = symbols) {
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
            output_string = output.join(', ')

            all_phrases.push(`<p> ${j+1}: ${input_string} \u{2192} ${output_string} </p>`)


          }

          return `<hr> <center> <p> Tested Compounds: <\p> ${all_phrases.join(' ')} </center>`
        }
        else{
          return "<hr> <center> No compounds have been tested yet </center>"
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
    output_string = output.join(', ')

    var tested = {
      type: jsPsychHtmlButtonResponse,
      stimulus: '<p style="font-size:48px; color:red;">GREEN</p>',
      choices: ['Red', 'Green', 'Blue'],
      prompt: "<p>What color is the ink?</p>"
    }

    return `<p> ${input_string} \u{2192} ${output_string} </p>`;

  },
  choices: ["Continue", "I'm done"]

}


var loop_node = {
    timeline: [trial, output_trial],
    loop_function: function(data, c = count){
      //console.log(c)
        if(c > 9 || data["trials"][1]["response"] == 1){
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
          output_string = output.join(', ')

          all_phrases.push(`<p> ${j+1}: ${input_string} \u{2192} ${output_string} </p>`)


        }

        return `<p> Congratulations! <br> The experimenter will now hand you a piece of paper to record what you think each compound produces. <br>
        For reference, below are the compounds you tested. </p> <hr> 
        <center> <p> Tested Compounds: <\p> ${all_phrases.join(' ')} </center>`
      },

  //"<p> Congratulations! The experimenter will instruct you on what to do next</p>",
  choices: ["Complete"]
};
timeline.push(end);



/* start the experiment */
//jsPsych.run([timeline, loop_node]);
jsPsych.run(timeline);
