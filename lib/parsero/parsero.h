// Copyright @mariokostelac
#include <cstdlib>
#include <cassert>
#include <unistd.h>
#include <vector>
#include <queue>
#include <string>
#include <sstream>

namespace parsero {

    struct option_t {
        std::string format;
        std::string description;
        void (*callback)(char *);

        option_t(std::string format, std::string description, void (*callback)(char *)) :
            format(format), description(description), callback(callback) {}
    };

    struct argument_t {
        std::string name;
        void (*callback)(char *);
        bool list;

        argument_t(std::string name, void (*callback)(char *), bool list)
            : name(name), callback(callback), list(list) {}
    };

    std::vector<option_t> options;
    std::vector<argument_t> arguments;
    std::string header, footer;

    void help(char *name) {

        fprintf(stderr, "SYNOPSIS\n");

        // print call line
        fprintf(stderr, "\t%s [options] ", name);
        for (auto curr_arg : arguments) {
            fprintf(stderr, "<%s%s> ", curr_arg.name.c_str(), curr_arg.list == true ? "..." : "");
        }
        fprintf(stderr, "\n\n");

        fprintf(stderr, "DESCRIPTION\n");
        if (header.length() > 0) {
            fprintf(stderr, "%s\n\n", header.c_str());
        }

        // print options
        if (options.size() > 0) {
            fprintf(stderr, "The following options are available:\n");

            for (auto curr_option : options) {
                fprintf(stderr, "\t-%c\t %s\n", curr_option.format[0], curr_option.description.c_str());
            }
        }

        if (footer.length() > 0) {
            fprintf(stderr, "\n%s\n", footer.c_str());
        }

        fprintf(stderr, "\n");
    }

    void set_header(std::string h) {
        header = h;
    }

    void set_footer(std::string f) {
        footer = f;
    }

    void add_option(std::string format, std::string description, void (*callback)(char *)) {

        assert(format.length() == 1 || (format.length() == 2 && format[1] == ':'));
        assert(callback != NULL);

        options.push_back(option_t(format, description, callback));
    }

    void add_argument(std::string name, void (*callback)(char *)) {

        assert(callback != NULL);

        arguments.push_back(argument_t(name, callback, false));
    }

    void add_arguments_list(std::string name, void (*callback)(char *)) {

        assert(callback != NULL);

        arguments.push_back(argument_t(name, callback, true));
    }

    int parse(int argc, char **argv) {

        // prepare options string ro getopt
        std::ostringstream format_stream;

        for (auto option : options) {
            format_stream << option.format;
        }

        std::string options_format = format_stream.str();

        // process options
        int o;
        while ((o = getopt(argc, argv, options_format.c_str())) != -1) {

            if (o == 'h') {
                help(argv[0]);
                exit(0);
            }

            // yep, i know it is not the fastest implementation, but we don't have 100+ options
            for (auto option : options) {
                if (option.format[0] == o) {
                    option.callback(optarg);
                }
            }
        }

        // break if we have more argument definitions than arguments
        if ((int) arguments.size() > argc - optind) {
            help(argv[0]);
            exit(1);
        }

        if (arguments.size() == 0) return optind;

        // process arguments
        size_t arg_i = 0;
        while (arg_i < arguments.size() && optind < argc) {
            auto curr_arg = arguments[arg_i];
            curr_arg.callback(argv[optind++]);

            if (!curr_arg.list) arg_i++;
        }

        return optind;
    }
};
