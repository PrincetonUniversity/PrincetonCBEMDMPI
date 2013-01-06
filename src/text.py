# This module contains general functions that have to deal with text
# manipulation.

def answer_question(question, type, print_answers = [], permit_answers = []):
    """This function is designed to be used in interactive programs to answer
    questions.  It asks the question and ensures that an acceptable answer is
    received."""
    accept = "no"
    while accept == "no":
        # Get an answer to the question
        print question
        answer = raw_input()
        flag = 0
        # If it is supposed to be an integer, ensure that it is
        if type == "integer":
            line = "You did not answer with an integer.  Please do so. "
            try:
                a1 = int(answer)
                a2 = float(answer)
                if a1 != a2:
                    print line
                    flag = 1
                else:
                    answer = a1
            except ValueError:
                print line
                flag = 1
        # If it is supposed to be a floating point number, ensure that it is
        elif type == "floating":
            try:
                answer = float(answer)
            except ValueError:
                print "You did not answer with a number.  Please do so. "
                flag = 1
        else:
            pass

        if flag == 0:
            # If no desired answers were provided in the function call, then the
            # answer is automatically acceptable
            if print_answers == [] and permit_answers == []:
                accept = "yes"
            # If the answer is in those provided in the function call, keep it
            elif answer in print_answers or answer in permit_answers:
                accept = "yes"
            # Otherwise provide information to the user about what they did
            # wrong
            else:
                print "You did not answer correctly.  Please try again."
                line = "Possible answers include: "
                if len(print_answers) == 1:
                    if print_answers[0] == ' ':
                        line += "' '"
                    else:
                        line += str(print_answers[0])
                    print line
                elif len(print_answers) == 2:
                    if print_answers[0] == ' ':
                        line += "' ' and "
                    else:
                        line += str(print_answers[0]) + " and " 
                    if print_answers[1] == ' ':
                        line += "' '"
                    else:
                        line += str(print_answers[1])
                elif len(print_answers) > 2:
                    last = print_answers[-1]
                    for answer in print_answers:
                        if answer == ' ':
                            use = "' '"
                        else:
                            use = str(answer)
                        if answer != last:
                            line += use + ", "
                        else:
                            line += "and " + use
                    print line
                else:
                    pass
    return answer
