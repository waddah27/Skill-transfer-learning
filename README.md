# Skill-transfer-learning
estimation of muscle generated forces from EMG datasets (hand gestures records)
Here I am using Hill-type model to represent human upper extrimity muscles activity which where recorded by EMG(Electromyugraphy) electodes alongside with the associated muscle forces under a set of hand gestures.
The goal is to develop a ML Neuromusculoskeletal (NMS) model that estimates the generated muscle forces from recorded EMG data. (I used the dataset recorded from nature.com to train the designed model)
The model is subject related (Parametric mathimatical method) depends on Hill-type representation of muscle contraction.
the MAIN.m file is the main notebook of the project.
the NMS.m is the Neuromusculoskeletal model
DataProc is the function for EMG signal preprocessing
notchfilter to filter a band of unwanted frequencies in the signal (noise).
