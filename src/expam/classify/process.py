from expam.process.manager import ControlCenter


class ExpamClassifierProcesses(ControlCenter):
    @classmethod
    def from_method_dict(cls, logging_dir, config):
        base_center = super(ExpamClassifierProcesses, cls).from_dict(logging_dir, config)

        base_arguments = {
            "group_name": "group_name",
            "workers": "workers",
            "child_statuses": "_child_statuses",
            "phases": "phases",
            "phase_queues": "phase_queues",
            "child_queues": "child_queues",
            "children": "_children",
            "processors": "_processors",
            "transitions": "_transitions",
            "timeout": "_timeout",
        }
        base_attributes = {
            attr: getattr(base_center, attr_reference)
            for attr, attr_reference in base_arguments.items()
        }

        # Create child class from this instance.
        control_center = ExpamClassifierProcesses(logging_dir=logging_dir, **base_attributes)

        control_center.set_methods(
            {
                "classify": {
                    "processor": None,
                    "transition": None
                }
            }
        )

        return control_center

